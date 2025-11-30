cmake_minimum_required(VERSION 3.29)
project(tlb)

set(CMAKE_CXX_STANDARD 20)

add_executable(tlb main.cpp)

target_compile_options(tlb PRIVATE "-O0")
