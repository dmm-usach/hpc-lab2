#!/usr/bin/env bash
# Ejecuta una corrida. Uso: ./scripts/run.sh P T [N] [pasos] [cada_cuanto_vtk] [prefijo]
set -euo pipefail
cd "$(dirname "$0")/.."

P="${1:-4}"
T="${2:-4}"
N="${3:-256}"
STEPS="${4:-1000}"
OUT="${5:-0}"
PREFIX="${6:-gas}"

if [[ ! -x build/gas ]]; then
    echo "Falta build/gas. Corré ./scripts/build.sh primero." >&2
    exit 1
fi

OMP_NUM_THREADS="$T" mpirun -np "$P" --bind-to none \
    ./build/gas "$N" "$STEPS" "$T" "$OUT" "$PREFIX"
