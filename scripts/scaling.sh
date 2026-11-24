#!/usr/bin/env bash
# Estudio de scaling. Uso: N=256 STEPS=200 ./scripts/scaling.sh
# Variables: P_VALUES, T_VALUES, N, STEPS, OUT
set -euo pipefail
cd "$(dirname "$0")/.."

N="${N:-256}"
STEPS="${STEPS:-200}"
OUT="${OUT:--1}"
P_VALUES="${P_VALUES:-1 2 4 8}"
T_VALUES="${T_VALUES:-1 2 4 8}"

mkdir -p results
CSV="results/scaling_$(date +%Y%m%d_%H%M%S).csv"
echo "P,T,threads,time_loop_s,compute_max_s,comm_max_s,mass0,mass1" > "$CSV"

for P in $P_VALUES; do
    for T in $T_VALUES; do
        total=$((P * T))
        echo "== P=$P T=$T (threads totales: $total) =="
        out=$(OMP_NUM_THREADS="$T" mpirun -np "$P" --bind-to none \
              ./build/gas "$N" "$STEPS" "$T" "$OUT" 2>&1)
        echo "$out"
        tl=$(sed -n 's/.*time_loop=\([0-9.]*\).*/\1/p' <<<"$out")
        cm=$(sed -n 's/.*compute_max=\([0-9.]*\).*/\1/p' <<<"$out")
        cc=$(sed -n 's/.*comm_max=\([0-9.]*\).*/\1/p' <<<"$out")
        m0=$(sed -n 's/.*mass0=\([0-9]*\).*/\1/p' <<<"$out")
        m1=$(sed -n 's/.*mass1=\([0-9]*\).*/\1/p' <<<"$out")
        echo "$P,$T,$total,$tl,$cm,$cc,$m0,$m1" >> "$CSV"
    done
done

echo "CSV: $CSV"
