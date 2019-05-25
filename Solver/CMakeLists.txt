include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/Domain/include)

add_library(Solver ${CMAKE_CURRENT_SOURCE_DIR}/src/Solver.cpp)
target_link_libraries(Solver Domain)

add_subdirectory(unit_test)