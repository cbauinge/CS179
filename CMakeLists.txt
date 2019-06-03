cmake_minimum_required(VERSION 3.10)
project(CS179 CXX)

#Set Standard
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CUDA)
if (CUDA_FOUND)
    add_definitions(-DUSE_CUDA)
    enable_language(CUDA)
    include_directories(SYSTEM ${CUDA_INCLUDE_DIRS} ${CUDA_TOOLKIT_INCLUDE} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(GTEST_PATH ${CMAKE_SOURCE_DIR}/../googletest)
add_subdirectory(${GTEST_PATH} ${CMAKE_BINARY_DIR}/gtest)

set(EIGEN_PATH ${CMAKE_SOURCE_DIR}/../eigen)
add_subdirectory(${EIGEN_PATH} ${CMAKE_BINARY_DIR}/eigen)

include_directories(SYSTEM ${GTEST_PATH}/googletest/include)
include_directories(SYSTEM ${EIGEN_PATH})
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

add_executable(Test Testframe.cpp)
target_link_libraries(Test Domain Reader Solver Writer)
add_custom_target(IntegrationTest Test)