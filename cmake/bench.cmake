file(GLOB BENCH_SOURCE bench/*.cpp)

option(BENCHMARK_BIG_O "Calculate Big O in benchmark" OFF)
option(BENCHMARK_100M "Run against 100M points" OFF)
option(BENCHMARK_10M "Run against 100M points" OFF)
set(BENCH_SOURCES ${BENCH_SOURCES} include/delaunator.cpp)
add_executable(bench-tests ${BENCH_SOURCES})
if(BENCHMARK_BIG_O)
    message("-- BENCHMARK_BIG_O=1")
    target_compile_definitions(bench-tests PUBLIC BENCHMARK_BIG_O=1)
endif()
if(BENCHMARK_100M)
    message("-- BENCHMARK_100M=1")
    target_compile_definitions(bench-tests PUBLIC BENCHMARK_100M=1)
endif()
if(BENCHMARK_10M)
    message("-- BENCHMARK_10M=1")
    target_compile_definitions(bench-tests PUBLIC BENCHMARK_10M=1)
endif()
