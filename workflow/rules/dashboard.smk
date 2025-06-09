from pathlib import Path

TREE = config["TREE"]
GIVEN_TAXONIUM = config.get("TAXONIUM_FILE", '')
if GIVEN_TAXONIUM:
    TAXONIUM_FILENAME = Path(GIVEN_TAXONIUM).name 
else:
    TAXONIUM_FILENAME = Path(TREE).stem + ".jsonl.gz"

rule process_taxonium:
    input:
        "build/wepp"
    output:
        jsonl=f"results/{{DIR}}/{TAXONIUM_FILENAME}"
    conda:
        "../envs/wepp.yml"
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"), 
        taxonium_jsonl_file=GIVEN_TAXONIUM,
        tree=TREE,
    shell:
        """ 
        if [ "{params.dashboard}" = "True" ]; then
            if [ "{params.taxonium_jsonl_file}" = '' ]; then
                echo "convert MAT file : data/{wildcards.DIR}/{params.tree} to {output.jsonl}"
                usher_to_taxonium --input data/{wildcards.DIR}/{params.tree} \
                    --output {output.jsonl} \
                    --name_internal_nodes -j src/Dashboard/taxonium_backend/config_public.json
            else
                cp data/{wildcards.DIR}/{params.taxonium_jsonl_file} {output.jsonl}
            fi
        else
            echo "Dashboard disabled."
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
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        bam_file="results/{DIR}/{FILE_PREFIX}_haplotype_reads.bam"
    conda:
        "../envs/wepp.yml"
    shell:
        """
        echo {params.bam_file}
        if [ "{params.dashboard}" = "True" ]; then

            if [ -f "{params.bam_file}" ]; then
                # Split BAM files by read groups
                echo "Splitting BAM file by read groups..."
                workflow/scripts/split_bams.sh {params.bam_file} ./results/{wildcards.DIR} {workflow.cores}

                rm {params.bam_file}
                rm {params.bam_file}.bai
            else
                echo "Splitting by read-groups already done!"
            fi
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
        "../envs/wepp.yml"
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        taxonium_jsonl_file=GIVEN_TAXONIUM,
        log=lambda wildcards: f"intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_run_tmp.txt",
        ref=config["REF"],
    shell:
        """
        if [ "{params.dashboard}" = "False" ]; then
            echo "removing .jsonl file that created. {input.taxonium_jsonl}"
            rm {input.taxonium_jsonl}
        fi

        echo "copying the reference file and indexing..." | tee -a {params.log}
        cp ./data/{wildcards.DIR}/{params.ref} ./results/{wildcards.DIR}/{params.ref} 
        samtools faidx ./results/{wildcards.DIR}/{params.ref}

        # Appending this project to database.
        python ./src/Dashboard/taxonium_backend/projects.py {wildcards.DIR} {input.taxonium_jsonl} ./results/{wildcards.DIR}/{params.ref}

        if [ "{params.dashboard}" = "True" ]; then
                # Start the Node.js server in the background
            if ss -tuln | grep ':8080' > /dev/null; then
                echo "Port 8080 is already in use. Exiting."
            else
                echo "creating uploads directory" | tee -a {params.log}
                if [ ! -d "./results/uploads" ]; then
                    mkdir ./results/uploads
                    cp ./src/Dashboard/data/NC_045512v2.* ./results/uploads/.
                fi
                
                echo "Installing backend dependencies..." | tee -a {params.log}
                cd src/Dashboard/taxonium_backend/ && yarn install && cd ../../../.

                echo "Starting the Node.js server..." | tee -a {params.log}

                # Get file size in MB
                MAX_MEM=$(du -m {input.taxonium_jsonl} | cut -f1)

                # Check if file is gzipped
                if [[ "{input.taxonium_jsonl}" == *.gz ]]; then
                    MAX_MEM=$(( MAX_MEM * 10 ))
                fi

                # Ensure at least 2048 MB, and add 1024 MB buffer
                MAX_MEM=$(( MAX_MEM < 2048 ? 2048 : MAX_MEM + 2048 ))

                echo "Allocating $MAX_MEM MB for Node.js server based on input file size..." | tee -a {params.log}

                node --expose-gc --max-old-space-size=$MAX_MEM src/Dashboard/taxonium_backend/server.js --port 8080 --data_file {input.taxonium_jsonl} &

                # Wait until port 8080 is open
                until ss -tuln | grep ':8080' > /dev/null; do
                    echo "Waiting for server to open port 8080..."
                    sleep 1
                done
            fi

            if ss -tuln | grep -w '80' > /dev/null; then
                echo "Port 80 is already in use. Exiting." | tee -a {params.log}
            else
                echo "Starting dashboard..." | tee -a {params.log}
                cp -r src/Dashboard/dashboard/dist/* /usr/share/nginx/html

                # create force symlink for static absolute paths for nginx to serve
                mkdir -p /srv/wepp
                ln -sfn "$(realpath ./results)" /srv/wepp/results

                cp src/Dashboard/nginx/nginx-react.conf /etc/nginx/sites-available/default
                nginx &

                until ss -tuln | grep -w '80' > /dev/null; do
                    echo "Waiting for nginx to start on port 80..."
                    sleep 1
                done
                
                echo "Dashboard can be accessed on browser on http://localhost:80 " | tee -a {params.log}
            fi

        else
            echo "Dashboard is disabled. Skipping Dashboard initialization."
        fi
        
        cp {params.log} {output}
        """