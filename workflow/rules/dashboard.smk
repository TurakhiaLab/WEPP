from pathlib import Path

TREE = config["TREE"]
GIVEN_TAXONIUM = config.get("TAXONIUM_FILE", '')
if GIVEN_TAXONIUM:
    TAXONIUM_FILENAME = Path(GIVEN_TAXONIUM).name 
else:
    TAXONIUM_FILENAME = Path(TREE).stem + ".jsonl.gz"

rule process_taxonium:
    input:
        str(BASE_DIR / "build/wepp")
    output:
        jsonl=f"results/{{DIR}}/{TAXONIUM_FILENAME}"
    conda:
        "../envs/dashboard.yml"
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"), 
        taxonium_jsonl_file=GIVEN_TAXONIUM,
        tree=TREE,
        clade_list=config["CLADE_LIST"],
    shell:
        """ 
        if [ "{params.dashboard}" = "True" ]; then
            if [ "{params.taxonium_jsonl_file}" = '' ]; then
                echo "convert MAT file : data/{wildcards.DIR}/{params.tree} to {output.jsonl}"
                usher_to_taxonium --input data/{wildcards.DIR}/{params.tree} \
                    --output {output.jsonl} \
                    --clade_types {params.clade_list} \
                    --name_internal_nodes
            else
                cp data/{wildcards.DIR}/{params.taxonium_jsonl_file} {output.jsonl}
            fi
        else
            echo "Dashboard disabled. Not converting pb file to jsonl.gz file."
            touch {output.jsonl}
        fi
        """

rule process_dashboard:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_run_tmp.txt",
        jsonl=f"results/{{DIR}}/{TAXONIUM_FILENAME}",
    output:
        log=temp("results/{DIR}/{FILE_PREFIX}_split_bam_log.txt")
    params:
        split_script = str(BASE_DIR / "workflow/scripts/split_bams.sh"),
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        bam_file="results/{DIR}/{FILE_PREFIX}_haplotype_reads.bam",
        haplotype_bam_file="{FILE_PREFIX}_haplotypes.bam",
        ref=config["REF"]
    conda:
        str(BASE_DIR / "workflow/envs/dashboard.yml")
    shell:
        """
        if [ "{params.dashboard}" = "True" ]; then
            if [ -f "{params.bam_file}" ]; then
                # Split BAM files by read groups
                echo "Splitting BAM file by read groups..."
                bash {params.split_script} {params.bam_file} ./results/{wildcards.DIR}/bams {workflow.cores}

                mv ./results/{wildcards.DIR}/{params.haplotype_bam_file} ./results/{wildcards.DIR}/bams/{params.haplotype_bam_file}
                mv ./results/{wildcards.DIR}/{params.haplotype_bam_file}.bai ./results/{wildcards.DIR}/bams/{params.haplotype_bam_file}.bai

                rm {params.bam_file}
                rm {params.bam_file}.bai
            else
                echo "Splitting by read-groups already done!"
            fi

            cp ./data/{wildcards.DIR}/{params.ref} ./results/{wildcards.DIR}/{params.ref} 
            samtools faidx ./results/{wildcards.DIR}/{params.ref}

        fi
        touch {output.log}
        """

rule dashboard_serve:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_run_tmp.txt",
        "results/{DIR}/{FILE_PREFIX}_split_bam_log.txt",
        taxonium_jsonl=f"results/{{DIR}}/{TAXONIUM_FILENAME}",
    output:
        temp("results/{DIR}/{FILE_PREFIX}_run.txt"),
    conda:
        str(BASE_DIR / "workflow/envs/dashboard.yml")
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        taxonium_jsonl_file=GIVEN_TAXONIUM,
        log=lambda wildcards: f"intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_run_tmp.txt",
        ref=config["REF"],
        find_port_script = str(BASE_DIR / "workflow/scripts/find_free_port.sh"),
        projects_py = str(BASE_DIR / "src/Dashboard/taxonium_backend/projects.py"),
        backend_dir = str(BASE_DIR / "src/Dashboard/taxonium_backend"),
        server_js = str(BASE_DIR / "src/Dashboard/taxonium_backend/server.js"),
        config_json = str(BASE_DIR / "src/Dashboard/taxonium_backend/config_public.json"),
        dist_dir = str(BASE_DIR / "src/Dashboard/dashboard/dist"),
        nginx_template = str(BASE_DIR / "src/Dashboard/nginx/wepp-nginx.conf.template"),
        mime_types = str(BASE_DIR / "src/Dashboard/nginx/mime.types")
    shell:
        """
        set -euo pipefail

        if [ "{params.dashboard}" = "False" ]; then
            rm -f {input.taxonium_jsonl}
        fi

        if [ "{params.dashboard}" = "True" ]; then

            export WEPP_DASHBOARD_PATH="$(pwd)/runtime"
            mkdir -p "$WEPP_DASHBOARD_PATH"

            # Check if running inside Docker
            if [ -f /.dockerenv ]; then
                IN_DOCKER=True
                export BACKEND_PORT=8080
            else
                source {params.find_port_script}
                IN_DOCKER=False
                export BACKEND_PORT=$(find_free_port)
            fi
            

            # ──────────────────────────────────────────────
            # Free up port $BACKEND_PORT if occupied
            # ──────────────────────────────────────────────
            if [ -f "$WEPP_DASHBOARD_PATH/dashboard.pid" ]; then
                DASH_PID=$(cat "$WEPP_DASHBOARD_PATH/dashboard.pid")
                if kill -0 "$DASH_PID" 2>/dev/null; then
                    echo -e "Existing dashboard running with PID $DASH_PID. Killing it safely."
                    kill -9 "$DASH_PID" || true
                else
                    echo -e "dashboard.pid found, but process $DASH_PID not running. Cleaning up file."
                fi
                rm -f "$WEPP_DASHBOARD_PATH/dashboard.pid"
            fi

            read FILENAME MAX_MEM < <( python {params.projects_py} \
                {wildcards.DIR} {input.taxonium_jsonl} ./results/{wildcards.DIR}/{params.ref})

            if [ ! -d "./results/uploads" ]; then
                echo "creating uploads directory" | tee -a {params.log}
                mkdir -p ./results/uploads
            fi

            echo "Installing backend dependencies..." | tee -a {params.log}
            (cd {params.backend_dir} && yarn install)

            echo -e "Starting the Node.js server..." | tee -a {params.log}

            if [[ "${{FILENAME}}" == *.gz ]]; then
                MAX_MEM=$(( MAX_MEM * 10 ))
            fi
            MAX_MEM=$(( MAX_MEM + 2048 ))
            echo -e "Allocating $MAX_MEM MB for Node.js server ..." | tee -a {params.log}

            node --expose-gc --max-old-space-size=$MAX_MEM \
                {params.server_js} \
                --port $BACKEND_PORT \
                --data_file {input.taxonium_jsonl} \
                --config_json {params.config_json} & 
            BACKEND_PID=$!

            echo -e "$BACKEND_PID" > "$WEPP_DASHBOARD_PATH/dashboard.pid"

            # ──────────────────────────────────────────────
            # Wait until Node.js server opens port $BACKEND_PORT
            # ──────────────────────────────────────────────
            until lsof -i :$BACKEND_PORT >/dev/null 2>&1; do
                # echo "Waiting for Node.js server to start..."
                sleep 1
            done

            # ──────────────────────────────────────────────
            # Start Nginx if not already in use
            # ──────────────────────────────────────────────

            if [ -f "$WEPP_DASHBOARD_PATH/nginx.pid" ]; then
                NGINX_PID=$(cat "$WEPP_DASHBOARD_PATH/nginx.pid")
                if kill -0 "$NGINX_PID" 2>/dev/null; then
                    echo "Existing nginx running with PID $NGINX_PID.  Killing it safely." | tee -a {params.log}
                    kill -9 "$NGINX_PID" || true
                    kill -9 $((NGINX_PID + 1)) || true
                else
                    echo "nginx.pid found, but process $NGINX_PID not running. Cleaning up file."
                    rm -f "$WEPP_DASHBOARD_PATH/nginx.pid"
                fi
            fi
                echo "Starting dashboard..." | tee -a {params.log}

                mkdir -p "$WEPP_DASHBOARD_PATH/Dashboard" "$WEPP_DASHBOARD_PATH/results"

                cp -r {params.dist_dir} "$WEPP_DASHBOARD_PATH/Dashboard"
                rm -rf "$WEPP_DASHBOARD_PATH/results"
                ln -sfn "$(realpath ./results)" "$WEPP_DASHBOARD_PATH/results"

                # Check if running inside Docker
                if [ "$IN_DOCKER" = "True" ]; then
                    export USER_DIRECTIVE="user root";
                    export FRONTEND_PORT=80
    
                else
                    source {params.find_port_script}
                    export USER_DIRECTIVE="# not user root";
                    export FRONTEND_PORT=$(find_free_port)
                fi

                envsubst '${{USER_DIRECTIVE}} ${{WEPP_DASHBOARD_PATH}} ${{FRONTEND_PORT}} ${{BACKEND_PORT}} ${{PWD}}' \
                < {params.nginx_template} \
                > "$WEPP_DASHBOARD_PATH/wepp-nginx.conf"

                cp {params.mime_types} $WEPP_DASHBOARD_PATH
                nginx -c $WEPP_DASHBOARD_PATH/wepp-nginx.conf &

                until lsof -i :$FRONTEND_PORT >/dev/null 2>&1; do
                    echo -e "Waiting for nginx to start on port $FRONTEND_PORT..."
                    sleep 1
                done
                if lsof -i :$FRONTEND_PORT >/dev/null 2>&1; then
                    echo -e "\n\n\nWORKFLOW COMPLETED! DASHBOARD IS RUNNING AT http://localhost:$FRONTEND_PORT (OR YOUR FORWARDED HOST PORT).\n\n\n"
                else
                    echo -e "\n\n\nWORKFLOW COMPLETED, BUT DASHBOARD NOT DETECTED ON PORT $FRONTEND_PORT. \n\n\n"
                fi

        else
            rm -f {input.taxonium_jsonl}
            echo -e "\n\n\nWorkflow completed! To run the dashboard, set DASHBOARD_ENABLED=True and rerun the workflow with --forcerun dashboard_serve.\n\n\n"

        fi

        cp {params.log} {output}
        """