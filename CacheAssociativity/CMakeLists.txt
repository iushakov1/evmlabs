cmake_minimum_required(VERSION 3.29)
project(CacheAssociativity)

set(CMAKE_CXX_STANDARD 20)

add_executable(CacheAssociativity main.cpp)

target_compile_options(CacheAssociativity PRIVATE "-O0")
