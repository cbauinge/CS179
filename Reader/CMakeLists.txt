include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/Domain/include)

add_library(Reader ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader.cpp)
target_link_libraries(Reader Domain)

add_subdirectory(unit_test)