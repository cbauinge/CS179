add_executable(DomainTest Testframe.cpp)
target_link_libraries(DomainTest Domain gtest_main)

add_test(DomainTest ${CMAKE_CURRENT_BINARY_DIR}/DomainTest)