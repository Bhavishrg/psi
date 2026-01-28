add_library(GraphDB
    utils/circuit.cpp
    utils/types.cpp
    utils/helpers.cpp
    graphdb/sharing.cpp
    graphdb/rand_gen_pool.cpp
    graphdb/offline_evaluator.cpp
    graphdb/online_evaluator_load_balanced.cpp)

target_include_directories(GraphDB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(GraphDB PUBLIC Boost::system EMPTool NTL GMP OpenMP::OpenMP_CXX rt tbb)
