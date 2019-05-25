#include "gtest/gtest.h"

#include "Reader.h"



class ReaderTest : public ::testing::Test
{
    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Reader reader_;
};


TEST(ReaderTest, constructor)
{
    ReaderTest r;
    r.reader_.read("test.bmp");
}