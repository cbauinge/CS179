cmake_minimum_required(VERSION 3.10)
project(CS179 CXX)

#Set Standard
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GTEST_PATH ${CMAKE_SOURCE_DIR}/../googletest)
add_subdirectory(${GTEST_PATH} ${CMAKE_BINARY_DIR}/gtest)

set(EIGEN_PATH ${CMAKE_SOURCE_DIR}/../eigen)
add_subdirectory(${EIGEN_PATH} ${CMAKE_BINARY_DIR}/eigen)

include_directories(SYSTEM ${GTEST_PATH}/googletest/include)
include_directories(${EIGEN_PATH})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Domain/include)
include_directories(${CMAKE_SOURCE_DIR}/Reader/include)
include_directories(${CMAKE_SOURCE_DIR}/Writer/include)
include_directories(${CMAKE_SOURCE_DIR}/Solver/include)


enable_testing()

add_subdirectory(Domain)
add_subdirectory(Reader)
add_subdirectory(Solver)
add_subdirectory(Writer)

add_executable(CS179 main.cpp)
target_link_libraries(CS179 Domain Reader Solver Writer)