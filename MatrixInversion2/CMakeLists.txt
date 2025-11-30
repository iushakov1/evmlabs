cmake_minimum_required(VERSION 3.29)
project(MatrixInversion2)

set(CMAKE_CXX_STANDARD 20)

add_executable(MatrixInversion2 main.cpp matrix.cpp)
target_compile_options(MatrixInversion2 PRIVATE "-O3")

enable_testing()

find_package(GTest REQUIRED)
find_package(BLAS REQUIRED)

add_executable(Tests matrix.cpp tests.cpp)

target_link_libraries(Tests
        GTest::gtest_main
        pthread
        ${BLAS_LIBRARIES}
)

target_link_libraries(MatrixInversion2
        ${BLAS_LIBRARIES}
)

include(GoogleTest)
gtest_discover_tests(Tests)
