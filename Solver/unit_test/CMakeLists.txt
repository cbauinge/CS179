add_definitions(-DDATA_PATH=${CMAKE_CURRENT_SOURCE_DIR})

add_executable(SolverTest Testframe.cpp)
target_link_libraries(SolverTest Solver gtest_main)

add_test(SolverTest SolverTest)