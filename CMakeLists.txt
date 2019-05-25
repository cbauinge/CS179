cmake_minimum_required(VERSION 3.0)
project(CS179 CXX)

#Set Standard
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GTEST_PATH ${CMAKE_SOURCE_DIR}/../googletest)

add_subdirectory(${GTEST_PATH} ${CMAKE_BINARY_DIR}/gtest)

include_directories(SYSTEM ${GTEST_PATH}/googletest/include)
include_directories(${CMAKE_SOURCE_DIR})

enable_testing()

add_subdirectory(Domain)
add_subdirectory(Reader)

add_executable(CS179 main.cpp)
target_link_libraries(CS179 Domain Reader)