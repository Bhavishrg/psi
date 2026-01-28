set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarks)

set(benchbin "")

macro(add_benchmark source_name)
    add_executable(${source_name} ${source_name}.cpp utils.cpp)
    target_link_libraries(${source_name} Boost::system Boost::program_options nlohmann_json::nlohmann_json GraphDB Threads::Threads NTL GMP EMPTool)
    list(APPEND benchbin ${source_name})
endmacro()


add_benchmark(graphiti_init)
add_benchmark(grapharo_init)
add_benchmark(bfs)
add_benchmark(insertV)
add_benchmark(insertE)
add_benchmark(modify)
add_benchmark(delete)
add_benchmark(naiveInsert)
add_benchmark(naiveDelete)
add_benchmark(naiveModify)

add_custom_target(benchmarks)
add_dependencies(benchmarks ${benchbin})