from pathlib import Path

TREE = config["TREE"]
TAXONIUM_PATH = config.get("USHER_TAXONIUM_FILE_PATH", '')
if TAXONIUM_PATH:
    TAXONIUM_FILENAME = Path(TAXONIUM_PATH).name 
else:
    TAXONIUM_FILENAME = Path(TREE).stem + ".jsonl"

rule process_taxonium:
    input:
        "build/wepp",
    output:
        f"results/{{DIR}}/{TAXONIUM_FILENAME}"
    params:
        taxonium_jsonl_file=TAXONIUM_PATH,
        tree=TREE,
    shell:
        """
         if [ "{params.taxonium_jsonl_file}" = '' ]; then
                echo "convert MAT file : data/{wildcards.DIR}/{params.tree} to {output}"
                usher_to_taxonium --input data/{wildcards.DIR}/{params.tree} \
                    --output {output} \
                    --name_internal_nodes -j src/Dashboard/taxonium_backend/config_public.json
            else
                echo "Taxonium file provided. Not converting MAT to .jsonl file"
                cp {params.taxonium_jsonl_file} {output}
            fi
        """

rule dashboard_serve:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb",
        "intermediate/{DIR}/{FILE_PREFIX}_run_tmp.txt",
        taxonium_jsonl=f"results/{{DIR}}/{TAXONIUM_FILENAME}",
    output:
        "results/{DIR}/{FILE_PREFIX}_run.txt"
    conda:
        "../envs/wepp.yml"
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        log=lambda wildcards: f"intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_run_tmp.txt",
        ref=config["REF"]
    shell:
        """
        echo "copying the reference file and indexing..." | tee -a {params.log}
        cp /workspace/data/{wildcards.DIR}/{params.ref} /workspace/results/{wildcards.DIR}/{params.ref}
        samtools faidx /workspace/results/{wildcards.DIR}/{params.ref}

        # Appending this project to database.
        python /workspace/src/Dashboard/taxonium_backend/projects.py {wildcards.DIR} {input.taxonium_jsonl} /workspace/results/{wildcards.DIR}/{params.ref}

        echo "Dashboard is enabled. creating taxonium file. {params.dashboard}"
        if [ "{params.dashboard}" = "True" ]; then
                # Start the Node.js server in the background
            if ss -tuln | grep ':8080' > /dev/null; then
                echo "Port 8080 is already in use. Exiting."
            else
                echo "creating uploads directory" | tee -a {params.log}
                if [ ! -d "/workspace/results/uploads" ]; then
                    mkdir /workspace/results/uploads
                    cp /workspace/src/Dashboard/data/NC_045512v2.* /workspace/results/uploads/.
                fi
                
                echo "Installing dashboard dependencies..." | tee -a {params.log}
                cd src/Dashboard/taxonium_backend/ && yarn install && cd /workspace

                echo "Starting the Node.js server..." | tee -a {params.log}

                node --expose-gc src/Dashboard/taxonium_backend/server.js --port 8080 --data_file {input.taxonium_jsonl} --integrated &

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
                cp src/Dashboard/nginx/nginx-react.conf /etc/nginx/sites-available/default
                nginx &

                until ss -tuln | grep -w '80' > /dev/null; do
                    echo "Waiting for nginx to start on port 80..."
                    sleep 1
                done
                
                echo "Dashboard serving on port 80..." | tee -a {params.log}
            fi

        else
            echo "Dashboard is disabled. Skipping input file processing."
        fi
        
        cp {params.log} {output}
        """

