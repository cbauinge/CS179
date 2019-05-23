set(GTEST_PATH ${CMAKE_SOURCE_DIR}/../googletest)

add_subdirectory(${GTEST_PATH} ${CMAKE_BINARY_DIR}/gtest)

include_directories(SYSTEM ${GTEST_PATH}/googletest/include)

add_executable(unit-test Testframe.cpp)
target_link_libraries(unit-test Domain gtest_main)

add_test(DomainTest unit-test)