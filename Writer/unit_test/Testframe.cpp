#include "Writer.h"
#include "Domain.h"
#include "gtest/gtest.h"

class WriterTest : public ::testing::Test
{
public:
    virtual void TestBody(){}    
    virtual void SetUp(){}
    virtual void TearDown(){}

    Writer writer_;

};


TEST(WriterTest, Test1)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {true, false, false};
    v.push_back(r1);
    std::vector<bool> r2 = {false, true, false};
    v.push_back(r2);
    v.push_back(r1);

    Domain dom{v};
    std::vector<double> result{3, 1.0};

    WriterTest wtest;
    wtest.writer_.Write("testoutput.csv", dom, result);
}