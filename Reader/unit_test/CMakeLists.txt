add_definitions(-DDATA_PATH=${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ReaderTest Testframe.cpp)
target_link_libraries(ReaderTest Reader gtest_main)

add_test(ReaderTest ReaderTest)