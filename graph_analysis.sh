#!/bin/bash
set -e

# graph_analysis.sh
# Runs a set of benchmarks (multiple graph sizes) using run.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SH="$SCRIPT_DIR/run.sh"

if [ ! -x "$RUN_SH" ]; then
    echo "Error: run.sh not found or not executable at $RUN_SH"
    exit 1
fi

SIZES=("1000:10000" "10000:100000" "100000:1000000" "1000000:10000000")
SIZES_GRAPHITI=("1000:10000" "10000:100000" "100000:1000000")

run_benchmark() {
    local bench="$1"
    shift
    local sizes=("${SIZES[@]}")
    if [ "$1" = "--sizes" ]; then
        shift
        sizes=("$@")
    fi
    for s in "${sizes[@]}"; do
        verts=${s%%:*}
        edges=${s##*:}
        echo "Running $bench with --num-verts $verts --num-edges $edges"
        if [ "$bench" = "add_vertices" ] || [ "$bench" = "add_edges" ]; then
            $RUN_SH "$bench" --num-vert "$verts" --num-edge "$edges" --latency 0 --random-inputs 1 --num-clients 10
        elif [ "$bench" = "graphiti_init" ]; then
            $RUN_SH "$bench" --num-vert "$verts" --num-edge "$edges" --latency 0
        elif [ "$bench" = "bfs" ]; then
            $RUN_SH "$bench" --num-vert "$verts" --num-edge "$edges" --iterations 2 --latency 0 --random-inputs 1
        elif [ "$bench" = "naiveInsert" ]; then
            $RUN_SH "$bench" --num-verts "$verts" --num-edges "$edges" --latency 0
        elif [ "$bench" = "naiveDelete" ] || [ "$bench" = "naiveModify" ]; then
            $RUN_SH "$bench" --num-verts "$verts" --num-edges "$edges" --latency 0
        else
            $RUN_SH "$bench" --num-vert "$verts" --num-edge "$edges" --latency 0 --random-inputs 1 --num-clients 10
        fi
        echo "Finished $bench $verts/$edges"
        echo "----------------------------------------"
    done
}

# echo "Starting graph analysis benchmarks"

# run_benchmark graphiti_init --sizes "${SIZES_GRAPHITI[@]}"
run_benchmark grapharo_init

# Run grapharo_init with varying client counts
echo "Running grapharo_init with varying client counts..."
for clients in 2 5 15; do
    echo "Running grapharo_init with --num-clients $clients --num-vert 100000 --num-edge 1000000"
    $RUN_SH grapharo_init --num-vert 100000 --num-edge 1000000 --num-clients "$clients" --latency 0 --random-inputs 1
    echo "Finished grapharo_init with $clients clients"
    echo "----------------------------------------"
done

run_benchmark insertV
run_benchmark insertE
run_benchmark delete
run_benchmark modify
run_benchmark bfs


echo "All benchmarks queued/completed. Check Results/ for logs."
