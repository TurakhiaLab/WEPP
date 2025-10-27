#!/usr/bin/env bash
find_free_port() {
    local port
    for port in $(seq 8000 9000); do
        if ! lsof -i :"$port" >/dev/null 2>&1; then
            echo "$port"
            return 0
        fi
    done
    return 1
}
