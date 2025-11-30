cmake_minimum_required(VERSION 3.29)
project(cacheCalc)

set(CMAKE_CXX_STANDARD 20)

add_executable(cacheCalc main.cpp)

target_compile_options(cacheCalc PRIVATE "-O0")