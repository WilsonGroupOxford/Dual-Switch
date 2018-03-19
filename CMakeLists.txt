cmake_minimum_required(VERSION 3.7)
project(dual_switch)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp simulation.cpp simulation.h customIO.cpp customIO.h network.cpp network.h node.cpp node.h crd.cpp crd.h vectorManip.cpp vectorManip.h minimise.cpp minimise.h)
add_executable(dual_switch ${SOURCE_FILES})