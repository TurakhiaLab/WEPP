import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from snakemake.exceptions import WorkflowError


def dashboard_is_enabled():
    return config.get("DASHBOARD_ENABLED", False) in (True, "True")


def dashboard_runtime_path(filename):
    return Path.cwd() / "runtime" / filename


def dashboard_pid_is_running(pid):
    if not isinstance(pid, int) or pid <= 0:
        return False

    try:
        state = subprocess.check_output(
            ["ps", "-o", "stat=", "-p", str(pid)],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return False

    return bool(state) and not state.startswith("Z")


def remove_owned_dashboard_pid_file(filename, expected_pid):
    pid_file = dashboard_runtime_path(filename)
    try:
        recorded_pid = pid_file.read_text().strip()
    except (FileNotFoundError, OSError):
        return

    if recorded_pid == str(expected_pid):
        try:
            pid_file.unlink()
        except FileNotFoundError:
            pass


def cleanup_dashboard_processes(backend_pid, nginx_pid):
    for pid in (nginx_pid, backend_pid):
        if dashboard_pid_is_running(pid):
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        if not any(dashboard_pid_is_running(pid) for pid in (backend_pid, nginx_pid)):
            break
        time.sleep(0.1)

    for pid in (nginx_pid, backend_pid):
        if dashboard_pid_is_running(pid):
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

    remove_owned_dashboard_pid_file("dashboard.pid", backend_pid)
    remove_owned_dashboard_pid_file("nginx.pid", nginx_pid)
    try:
        dashboard_runtime_path("dashboard.ready").unlink()
    except FileNotFoundError:
        pass


def monitor_dashboard():
    ready_file = dashboard_runtime_path("dashboard.ready")
    if not ready_file.is_file():
        return

    try:
        frontend_port, backend_pid, nginx_pid = ready_file.read_text().split()
        backend_pid = int(backend_pid)
        nginx_pid = int(nginx_pid)
        frontend_port = int(frontend_port)
    except (OSError, TypeError, ValueError):
        try:
            ready_file.unlink()
        except FileNotFoundError:
            pass
        raise WorkflowError("Dashboard lifecycle state is invalid.")

    if not dashboard_pid_is_running(backend_pid) or not dashboard_pid_is_running(nginx_pid):
        cleanup_dashboard_processes(backend_pid, nginx_pid)
        raise WorkflowError("Dashboard services exited before lifecycle handoff.")

    print(
        "\n\n\nDASHBOARD IS RUNNING AT "
        f"http://localhost:{frontend_port} (OR YOUR FORWARDED HOST PORT).\n"
        "Press Ctrl-C to stop the dashboard.\n\n\n",
        flush=True,
    )

    shutdown_requested = False

    def request_shutdown(signum, frame):
        nonlocal shutdown_requested
        shutdown_requested = True

    previous_sigint_handler = signal.signal(signal.SIGINT, request_shutdown)
    try:
        while not shutdown_requested:
            backend_running = dashboard_pid_is_running(backend_pid)
            nginx_running = dashboard_pid_is_running(nginx_pid)
            if shutdown_requested or not (backend_running and nginx_running):
                break
            time.sleep(0.5)
    finally:
        signal.signal(signal.SIGINT, previous_sigint_handler)

    if shutdown_requested:
        cleanup_dashboard_processes(backend_pid, nginx_pid)
        print("Dashboard stopped.", flush=True)
        return

    if not dashboard_pid_is_running(backend_pid):
        failed_service = "Dashboard backend"
    else:
        failed_service = "nginx"

    cleanup_dashboard_processes(backend_pid, nginx_pid)
    raise WorkflowError(f"{failed_service} exited unexpectedly.")


if dashboard_is_enabled() and "--no-hooks" in sys.argv:
    raise WorkflowError("The dashboard requires Snakemake hooks; remove --no-hooks.")


onstart:
    if dashboard_is_enabled():
        try:
            dashboard_runtime_path("dashboard.ready").unlink()
        except FileNotFoundError:
            pass


onsuccess:
    if dashboard_is_enabled():
        monitor_dashboard()

TREE = config["TREE"]
GIVEN_TAXONIUM = config.get("TAXONIUM_FILE") or ""
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
        sanitize_mat_newick=str(BASE_DIR / "workflow/scripts/sanitize_mat_newick.py"),
    shell:
        """ 
        if [ "{params.dashboard}" = "True" ]; then
            if [ "{params.taxonium_jsonl_file}" = '' ]; then
                echo "convert MAT file : data/{wildcards.DIR}/{params.tree} to {output.jsonl}"
                sanitized_mat=$(mktemp "${{TMPDIR:-/tmp}}/wepp_taxonium.XXXXXX.pb.gz")
                trap 'rm -f "$sanitized_mat"' EXIT
                python {params.sanitize_mat_newick} \
                    --input data/{wildcards.DIR}/{params.tree} \
                    --output "$sanitized_mat"
                usher_to_taxonium --input "$sanitized_mat" \
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
            rm -f "$WEPP_DASHBOARD_PATH/dashboard.ready"
            READY_FILE_TMP="$WEPP_DASHBOARD_PATH/dashboard.ready.$$"

            BACKEND_PID=""
            NGINX_PID=""

            pid_is_running() {{
                PID_STATE=$(ps -o stat= -p "$1" 2>/dev/null | tr -d '[:space:]')
                case "$PID_STATE" in
                    ""|Z*) return 1 ;;
                    *) return 0 ;;
                esac
            }}

            remove_owned_pid_file() {{
                PID_FILE=$1
                EXPECTED_PID=$2
                if [ -n "$EXPECTED_PID" ] && [ -f "$PID_FILE" ] && \
                    [ "$(cat "$PID_FILE")" = "$EXPECTED_PID" ]; then
                    rm -f "$PID_FILE"
                fi
            }}

            cleanup_dashboard() {{
                EXIT_CODE=$?
                trap - EXIT INT TERM

                for PID in "$NGINX_PID" "$BACKEND_PID"; do
                    if [ -n "$PID" ] && pid_is_running "$PID"; then
                        kill -TERM "$PID" 2>/dev/null || true
                    fi
                done

                for _ in 1 2 3 4 5; do
                    BACKEND_RUNNING=False
                    NGINX_RUNNING=False
                    if [ -n "$BACKEND_PID" ] && pid_is_running "$BACKEND_PID"; then
                        BACKEND_RUNNING=True
                    fi
                    if [ -n "$NGINX_PID" ] && pid_is_running "$NGINX_PID"; then
                        NGINX_RUNNING=True
                    fi
                    if [ "$BACKEND_RUNNING" = "False" ] && [ "$NGINX_RUNNING" = "False" ]; then
                        break
                    fi
                    sleep 1
                done

                for PID in "$NGINX_PID" "$BACKEND_PID"; do
                    if [ -n "$PID" ] && pid_is_running "$PID"; then
                        kill -KILL "$PID" 2>/dev/null || true
                    fi
                    if [ -n "$PID" ]; then
                        wait "$PID" 2>/dev/null || true
                    fi
                done

                remove_owned_pid_file "$WEPP_DASHBOARD_PATH/dashboard.pid" "$BACKEND_PID"
                remove_owned_pid_file "$WEPP_DASHBOARD_PATH/nginx.pid" "$NGINX_PID"
                rm -f "$WEPP_DASHBOARD_PATH/dashboard.ready" "$READY_FILE_TMP"
                exit "$EXIT_CODE"
            }}

            trap cleanup_dashboard EXIT
            trap 'exit 130' INT
            trap 'exit 143' TERM

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
                    kill -TERM "$DASH_PID" || true
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
                if ! kill -0 "$BACKEND_PID" 2>/dev/null; then
                    echo "Dashboard backend exited before opening port $BACKEND_PORT." | tee -a {params.log}
                    wait "$BACKEND_PID" || true
                    exit 1
                fi
                sleep 1
            done

            # ──────────────────────────────────────────────
            # Start Nginx if not already in use
            # ──────────────────────────────────────────────

            if [ -f "$WEPP_DASHBOARD_PATH/nginx.pid" ]; then
                NGINX_PID=$(cat "$WEPP_DASHBOARD_PATH/nginx.pid")
                if kill -0 "$NGINX_PID" 2>/dev/null; then
                    echo "Existing nginx running with PID $NGINX_PID.  Killing it safely." | tee -a {params.log}
                    kill -TERM "$NGINX_PID" || true
                else
                    echo "nginx.pid found, but process $NGINX_PID not running. Cleaning up file."
                    rm -f "$WEPP_DASHBOARD_PATH/nginx.pid"
                fi
            fi
                echo "Starting dashboard..." | tee -a {params.log}

                mkdir -p "$WEPP_DASHBOARD_PATH/Dashboard"

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
                nginx -g 'daemon off;' -c "$WEPP_DASHBOARD_PATH/wepp-nginx.conf" &
                NGINX_PID=$!

                until lsof -i :$FRONTEND_PORT >/dev/null 2>&1; do
                    if ! kill -0 "$NGINX_PID" 2>/dev/null; then
                        echo "nginx exited before opening port $FRONTEND_PORT." | tee -a {params.log}
                        wait "$NGINX_PID" || true
                        exit 1
                    fi
                    echo -e "Waiting for nginx to start on port $FRONTEND_PORT..."
                    sleep 1
                done
                if lsof -i :$FRONTEND_PORT >/dev/null 2>&1; then
                    cp {params.log} {output}
                    printf '%s %s %s\n' "$FRONTEND_PORT" "$BACKEND_PID" "$NGINX_PID" > "$READY_FILE_TMP"
                    mv "$READY_FILE_TMP" "$WEPP_DASHBOARD_PATH/dashboard.ready"
                else
                    echo -e "\n\n\nDASHBOARD NOT DETECTED ON PORT $FRONTEND_PORT. \n\n\n"
                    exit 1
                fi

        else
            rm -f {input.taxonium_jsonl}
            echo -e "\n\n\nWorkflow completed! To run the dashboard, set DASHBOARD_ENABLED=True and rerun the workflow with --forcerun dashboard_serve.\n\n\n"

        fi

        if [ "{params.dashboard}" = "True" ]; then
            trap - EXIT INT TERM
        else
            cp {params.log} {output}
        fi
        """
