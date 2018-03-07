cmake_minimum_required(VERSION 3.8)
project(dual_switch)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp simulation.cpp simulation.h customIO.cpp customIO.h)
add_executable(dual_switch ${SOURCE_FILES})