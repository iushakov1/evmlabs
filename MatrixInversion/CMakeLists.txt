cmake_minimum_required(VERSION 3.29)
project(MatrixInversion)

set(CMAKE_CXX_STANDARD 20)
add_executable(MatrixInversion main.cpp
        matrix.cpp)
target_compile_options(MatrixInversion PRIVATE "-O3" "-ftree-vectorize" "-fopt-info-vec-optimized")
#target_compile_options(MatrixInversion PRIVATE "-O3")
enable_testing()

find_package(GTest REQUIRED)

add_executable(Tests matrix.cpp
        tests.cpp)
target_link_libraries(Tests GTest::gtest_main pthread)

include(GoogleTest)
gtest_discover_tests(Tests)
target_compile_options(Tests PRIVATE "-O3")
