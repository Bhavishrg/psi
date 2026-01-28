# Grapharo

This directory contains the implementation of the Grapharo framework.
The protocol is implemented in C++17 and [CMake](https://cmake.org/) is used as the build system.

## External Dependencies
The following libraries need to be installed separately and should be available to the build system and compiler.

- [GMP](https://gmplib.org/)
- [NTL](https://www.shoup.net/ntl/) (11.0.0 or later)
- [Boost](https://www.boost.org/) (1.72.0 or later)
- [Nlohmann JSON](https://github.com/nlohmann/json)
- [EMP Tool](https://github.com/emp-toolkit/emp-tool)

### Docker
All required dependencies to compile and run the project are available through the docker image.
To build and run the docker image, execute the following commands from the root directory of the repository: 

```sh
# Build the grapharo Docker image.
#
# Building the Docker image requires at least 4GB RAM. This needs to be set 
# explicitly in case of Windows and MacOS.

docker build -t grapharo .

# Create and run a container.
#
# This should start the shell from within the container.
docker run -it -v $PWD:/code grapharo

# The following command changes the working directory to the one containing the 
# source code and should be run on the shell started using the previous command.
cd /code
```

## Compilation
The project uses [CMake](https://cmake.org/) for building the source code. 
To compile, run the following commands from the root directory of the repository:

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..

# The two main targets are 'benchmarks' and 'tests' corresponding to
# binaries used to run benchmarks and unit tests respectively.
make 
```


### Running Benchmarks

## Scripts

Two helper scripts are provided to run benchmarks and collect logs: [run.sh](run.sh) and [graph_analysis.sh](graph_analysis.sh).

### run.sh

Usage: `./../run.sh <benchmark_name> [benchmark_options...]`

**Available benchmarks:**

- `bfs` - Breadth-First Search benchmark
- `grapharo_init` - Grapharo initialization benchmark
- `graphiti_init` - Graphiti initialization benchmark
- `insertE` - Edge insertion benchmark
- `insertV` - Vertex insertion benchmark
- `modify` - Graph modification benchmark
- `delete` - Graph deletion benchmark

**Example usage:**

```sh
# From build directory, run BFS benchmark with 1000 vertices
./../run.sh bfs --num-verts 1000 -n 2

# Run insertE benchmark with custom parameters
./../run.sh insertE --num-verts 10000 --num-edges 50000 -n 3

# Run grapharo_init benchmark
./../run.sh grapharo_init --num-verts 100000 -n 2
```

### graph_analysis.sh

Simple wrapper to run all benchmarks reported in the paper and save logs.

**Usage:**

```sh
# From build directory:
./../graph_analysis.sh
```

## Results / Logs Layout

Both scripts store logs under the `Results/` directory with the following layout:

```
Results/<benchmark_name>/<num_verts>/<num_edges>/<num_clients>/
    party_0.log   # trusted party
    party_1.log
    party_2.log
    ...
    aggregate_stat.log
```

The `aggregate_stat.log` file reports the aggregated runtime and communication statistics.

The included Python scripts under `pythonScripts/getAggStat.py` can be used to aggregate and convert these logs into human-friendly tables.
