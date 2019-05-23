#include "Domain.h"
#include "gtest/gtest.h"


class DomainTest : public ::testing::Test
{
protected:
    virtual void SetUp(){}
    virtual void TearDown(){}

    Domain domain_;
};

TEST(DomainTest, Constructor)
{
    std::vector<std::vector<bool> > v{true};
    //Domain{v};
}