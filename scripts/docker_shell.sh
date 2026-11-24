#!/usr/bin/env bash
# Abre el contenedor con el repo montado en /work.
# Uso: ./scripts/docker_shell.sh bash -lc './scripts/build.sh && ./scripts/run.sh 4 4'
set -euo pipefail
cd "$(dirname "$0")/.."

TTY=()
if [[ -t 0 && -t 1 ]]; then
    TTY=(-it)
fi

docker run --rm "${TTY[@]}" \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e TMPDIR=/tmp \
    -v "$PWD":/work \
    -w /work \
    hpc-lab2:latest "$@"
