configfile: config.get("config_path", "config/config.yaml")

rule process_taxonium:
    input:
        "build/wepp",
    output:
        "results/{DIR}/"+config['TREE']+'.jsonl'
    params:
        taxonium_jsonl_file=config.get("USHER_TAXONIUM_FILE_PATH", ''),
        tree=config["TREE"],
    shell:
        """
         if [ "{params.taxonium_jsonl_file}" = '' ]; then
                
                echo "data/{wildcards.DIR}/{params.tree} {output}"
                usher_to_taxonium --input data/{wildcards.DIR}/{params.tree} \
                    --output {output} \
                    --name_internal_nodes -j src/Dashboard/taxonium_backend/config_public.json
            else
                echo "Dashboard is enabled and taxonium file provided."
                cp {params.taxonium_jsonl_file} {output}
            fi
        """

rule taxonium_backend:
    input:
        # "build/wepp",
        "results/{DIR}/"+config['TREE']+".jsonl",
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb",
        "intermediate/{DIR}/{FILE_PREFIX}_run_tmp.txt"
    output:
        "results/{DIR}/{FILE_PREFIX}_run.txt"
    params:
        dashboard=config.get("DASHBOARD_ENABLED", "false"),
        log=lambda wildcards: f"intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_run_tmp.txt",
        ref=config["REF"]
    shell:
        """
        echo "copying the reference file and indexing..." | tee -a {params.log}
        cp data/{wildcards.DIR}/{params.ref} ./results/{wildcards.DIR}/{params.ref}
        samtools faidx ./results/{wildcards.DIR}/{params.ref}


        # creating or appending projects.json
        python src/Dashboard/taxonium_backend/projects.py {wildcards.DIR} {input[0]} ./results/{wildcards.DIR}/{params.ref}


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

                node --expose-gc src/Dashboard/taxonium_backend/server.js --port 8080 --data_file {input[0]} --integrated &

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
        
        mv {params.log} {output}
        """