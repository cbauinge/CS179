include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Domain ${CMAKE_CURRENT_SOURCE_DIR}/src/Domain.cpp)

add_subdirectory(unit_test)