configfile: config.get("config_path", "config/config.yaml")

rule process_taxonium:
    input:
        "build/wbe",
    output:
         "intermediate/{dataset}/dashboard/taxonium.jsonl.gz"
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        taxonium_jsonl_file=config.get("USHER_TAXONIUM_FILE_PATH", ''),
        tree=config["TREE"],
    shell:
        """
        if [ "{params.dashboard}" = "True" ]; then
            if [ "{params.taxonium_jsonl_file}" = '' ]; then
                
                echo "Dashboard is enabled. creating taxonium file"

                usher_to_taxonium --input data/{wildcards.dataset}/{params.tree} \
                    --output {output} \
                    --name_internal_nodes -j src/dashboard/taxonium_backend/config_public.json
            else
                echo "Dashboard is enabled and taxonium file provided."
                cp {params.taxonium_jsonl_file} {output}
            fi
        else
            echo "Dashboard is disabled. Skipping input file processing."
        fi
        """

rule taxonium_backend:
    input:
        # "build/wbe",
        "intermediate/{dataset}/dashboard/taxonium.jsonl.gz",
        "intermediate/{dataset}/{file_prefix}_depth.tsv",
        "intermediate/{dataset}/{file_prefix}_corrected_variants.tsv",
        "intermediate/{dataset}/{file_prefix}_reads.pb",
        "intermediate/{dataset}/{file_prefix}_run_tmp.txt"
    output:
        "results/{dataset}/{file_prefix}_run.txt"
    params:
        log=lambda wildcards: f"intermediate/{wildcards.dataset}/{wildcards.file_prefix}_run_tmp.txt"
    shell:
        """
        # Start the Node.js server in the background
        if ss -tuln | grep ':8080' > /dev/null; then
            echo "Port 8080 is already in use. Exiting."
        else
            echo "copying the reference file and indexing..." | tee -a {params.log}
            cp data/{wildcards.dataset}/{config[REF]} ./results/{wildcards.dataset}/NC_045512v2.fa

            samtools faidx ./results/{wildcards.dataset}/NC_045512v2.fa

            echo "creating uploads directory" | tee -a {params.log}
            if [ ! -d "/workspace/results/uploads" ]; then
                mkdir /workspace/results/uploads
                cp /workspace/src/dashboard/data/NC_045512v2.* /workspace/results/uploads/.
            fi

            echo "Installing dashboard dependencies..." | tee -a {params.log}
            cd src/dashboard/taxonium_backend/ && yarn install && cd /workspace

            echo "Starting the Node.js server..." | tee -a {params.log}

            node src/dashboard/taxonium_backend/server.js --port 8080 --data_file {input[0]} --integrated &

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
            cp -r src/dashboard/dashboard/dist/* /usr/share/nginx/html
            cp src/dashboard/nginx/nginx-react.conf /etc/nginx/sites-available/default
            nginx &

            until ss -tuln | grep -w '80' > /dev/null; do
                echo "Waiting for nginx to start on port 80..."
                sleep 1
            done
            
            echo "Dashboard serving on port 80..." | tee -a {params.log}
        fi

        mv {params.log} {output}
        """