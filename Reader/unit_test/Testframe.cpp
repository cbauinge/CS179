#include "gtest/gtest.h"

#include "Reader.h"

#include <exception>
#include <iostream>
#include <string>


#ifndef DATA_PATH //this should be defined in the CMakeLists.txt in the same folder as this file.
    std::cout << "no data path defined" << std::endl;
    #define DATA_PATH "."
#endif

//weird construct necessary to get the value in a maro into a string
#define STR(s) STRTMP(s)
#define STRTMP(s) #s


class ReaderTest : public ::testing::Test
{
public:

    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Reader reader_;
};


TEST(ReaderTest, Test1)
{
    ReaderTest r;
    std::string path(STR(DATA_PATH));
    std::cout << "path = " << path << std::endl;
    std::string file1 = path + std::string("/test1.bmp");
    Domain d = r.reader_.Read(file1.c_str());
    EXPECT_EQ(0, d.GetPoints()[0][0]);
    EXPECT_EQ(1, d.GetPoints()[4][4]);
}