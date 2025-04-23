cmake_minimum_required (VERSION 3.16)
project (Benchmark_new_for_hybrid)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

find_package ( T8CODE 4.0 REQUIRED)

add_library(Benchmark_new_for_hybrid benchmark.cxx)

target_link_libraries (Benchmark_new_for_hybrid PRIVATE T8CODE::T8)

target_include_directories(Benchmark_new_for_hybrid
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

add_executable(benchmark benchmark.cxx)

target_link_libraries(benchmark PRIVATE T8CODE::T8)
