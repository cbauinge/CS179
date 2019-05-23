cmake_minimum_required(VERSION 3.0)
project(CS179 CXX)

#Set Standard
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)

add_executable(CS179 main.cpp)