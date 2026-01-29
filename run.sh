#!/bin/bash
# set -x

# Usage: ./run.sh <benchmark_name> [benchmark_options...]
# Example: ./run.sh grouppropagate -l 0.5 --t1-size 8 --t2-size 10
# Example: ./run.sh reconstruction -l 0.5 -i 10 --use-pking true
# Example: ./run.sh mult -l 100.0 -i 30

if [ $# -lt 1 ]; then
    echo "Usage: $0 <benchmark_name> [benchmark_options...]"
    echo ""
    echo "Available benchmarks (./benchmarks/):"
    ls -1 ./benchmarks/ 2>/dev/null | sed 's/^/  - /' || echo "  (directory not found)"
    echo ""
    echo "Available test primitives (./test_primitives/):"
    ls -1 ./test_primitives/ 2>/dev/null | sed 's/^/  - /' || echo "  (directory not found)"
    echo ""
    echo "Example: $0 grouppropagate -l 0.5 --t1-size 8 --t2-size 10"
    echo "Example: $0 reconstruction -l 0.5 -i 10 --use-pking true"
    exit 1
fi

# Get benchmark name from first argument
BENCHMARK_NAME="$1"
shift  # Remove first argument, leaving only the options

# Validate benchmark exists - check both benchmarks and test_primitives directories
BENCHMARK_PATH=""
if [ -f "./benchmarks/$BENCHMARK_NAME" ]; then
    BENCHMARK_PATH="./benchmarks/$BENCHMARK_NAME"
    BENCHMARK_TYPE="benchmarks"
elif [ -f "./test_primitives/$BENCHMARK_NAME" ]; then
    BENCHMARK_PATH="./test_primitives/$BENCHMARK_NAME"
    BENCHMARK_TYPE="test_primitives"
fi

if [ -z "$BENCHMARK_PATH" ]; then
    echo "Error: Benchmark '$BENCHMARK_NAME' not found in ./benchmarks/ or ./test_primitives/"
    echo ""
    echo "Available benchmarks in ./benchmarks/:"
    ls -1 ./benchmarks/ 2>/dev/null | sed 's/^/  - /' || echo "  (directory not found)"
    echo ""
    echo "Available test primitives in ./test_primitives/:"
    ls -1 ./test_primitives/ 2>/dev/null | sed 's/^/  - /' || echo "  (directory not found)"
    exit 1
fi

# Capture all remaining arguments as benchmark options (preserve spacing)
BENCHMARK_OPTS=("$@")

# Extract options we need for directory naming
players=2
num_verts="unspecified_verts"
num_edges="unspecified_edges"
num_clients=""

for ((i=0; i<${#BENCHMARK_OPTS[@]}; i++)); do
    case "${BENCHMARK_OPTS[$i]}" in
        -n|--num-parties)
            if (( i + 1 < ${#BENCHMARK_OPTS[@]} )); then
                players="${BENCHMARK_OPTS[$((i+1))]}"
            fi
            ;;
        --num-vert|--num-verts)
            if (( i + 1 < ${#BENCHMARK_OPTS[@]} )); then
                num_verts="${BENCHMARK_OPTS[$((i+1))]}"
            fi
            ;;
        --num-edge|--num-edges)
            if (( i + 1 < ${#BENCHMARK_OPTS[@]} )); then
                num_edges="${BENCHMARK_OPTS[$((i+1))]}"
            fi
            ;;
        --num-clients)
            if (( i + 1 < ${#BENCHMARK_OPTS[@]} )); then
                num_clients="${BENCHMARK_OPTS[$((i+1))]}"
            fi
            ;;
    esac
done

if [[ -z "$num_clients" ]]; then
    num_clients="$players"
fi

# If -n/--num-parties not found in options, add default
if [[ ! " ${BENCHMARK_OPTS[*]} " =~ " -n " ]] && [[ ! " ${BENCHMARK_OPTS[*]} " =~ " --num-parties " ]]; then
    BENCHMARK_OPTS+=("-n" "$players")
fi

# Create results directory structure: Results/<benchmark_type>/<benchmark_name>/<num_verts>/<num_edges>/<num_clients>
dir=$PWD/../Results/$BENCHMARK_TYPE/$BENCHMARK_NAME/$num_verts/$num_edges/$num_clients

# Clean up old results for this benchmark and party configuration
# rm -rf $dir

echo "Running benchmark: $BENCHMARK_NAME (from $BENCHMARK_TYPE)"
echo "Number of players: $players"
echo "Benchmark options: ${BENCHMARK_OPTS[*]}"
echo "Results directory: $dir"
echo ""

for rounds in 1
do
    for party in $(seq 1 $players)
    do
        logdir=$dir
        mkdir -p "$logdir"
        log=$logdir/party_$party.log
        tplog=$logdir/party_0.log

        # Run benchmark for each party with --localhost option
        "$BENCHMARK_PATH" "${BENCHMARK_OPTS[@]}" --localhost -p "$party" 2>&1 | cat > "$log" &
        codes[$party]=$!
    done
    
    # Run benchmark for party 0 (trusted party)
    "$BENCHMARK_PATH" "${BENCHMARK_OPTS[@]}" --localhost -p 0 2>&1 | cat > "$tplog" & 
    codes[0]=$!
    
    # Wait for all parties to complete
    for party in $(seq 0 $players)
    do
        wait ${codes[$party]} || return 1
    done
done

echo "Benchmark execution completed. Logs saved to: $logdir"
echo ""

# Run aggregation script if it exists
if [ -f "/code/pythonScripts/getAggStat.py" ]; then
    echo "Running aggregation script..."
    python3 /code/pythonScripts/getAggStat.py "$logdir/"
fi
