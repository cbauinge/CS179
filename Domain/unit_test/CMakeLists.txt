add_executable(DomainTest Testframe.cpp)
target_link_libraries(DomainTest Domain gtest_main)

add_test(test DomainTest)