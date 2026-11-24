#!/usr/bin/env bash
# Compila el proyecto (dentro del contenedor o en una máquina con MPI+BSPonMPI).
set -euo pipefail
cd "$(dirname "$0")/.."

cmake -S . -B build -DCMAKE_BUILD_TYPE="${BUILD_TYPE:-Release}"
cmake --build build -j"${JOBS:-$(nproc)}"
