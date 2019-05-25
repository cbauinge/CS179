include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/Domain/include)

add_library(Writer ${CMAKE_CURRENT_SOURCE_DIR}/src/Writer.cpp)
target_link_libraries(Writer Domain)

add_subdirectory(unit_test)