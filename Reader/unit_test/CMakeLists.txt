add_executable(ReaderTest Testframe.cpp)
target_link_libraries(ReaderTest Reader gtest_main)

add_test(test ReaderTest)