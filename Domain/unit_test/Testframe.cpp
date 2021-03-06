#include "Domain.h"
#include "Vec2D.h"
#include "BoundaryDataGenerator.h"
#include "gtest/gtest.h"
#include <cstdlib>

double testfunction(double a, double b)
{
    return 1.0;
}


class DomainTest : public ::testing::Test
{
public:
    DomainTest(std::vector<std::vector<bool> > points) : 
        domain_(points)
    {}

    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Domain domain_;
};


class Vec2DTest : public ::testing::Test
{
    Vec2DTest(double x, double y) : vec_{x, y} {}

    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Vec2D vec_;
};


class BoundaryDataGeneratorTest : public ::testing::Test
{
public:
    BoundaryDataGeneratorTest(double (*func)(double, double)) : generator_(func) {}

    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    BoundaryDataGenerator generator_;
};

class BoundaryTest : public ::testing::Test
{
public:
    BoundaryTest(std::vector<std::vector<bool> > points, double h) : boundary_(points, h) {}

    bool IsBoundary(int i, int j, std::vector<std::vector<bool> > points) {return boundary_.IsBoundary(i, j, points);}

    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Boundary boundary_;
};


TEST(Vec2DTest, Constructor)
{
    Vec2D v{1, 2};
    EXPECT_EQ(1, v.x());
    EXPECT_EQ(2, v.y());
}

TEST(Vec2DTest, Addition)
{
    Vec2D v1{1, 2};
    Vec2D v2{-2, 1};
    Vec2D v3 = v1 + v2;
    EXPECT_EQ(-1, v3.x());
    EXPECT_EQ(3, v3.y());
}

TEST(Vec2DTest, AssignmentAddition)
{
    Vec2D v1{1, 2};
    Vec2D v2{-2, 1};
    v1 += v2;
    EXPECT_EQ(-1, v1.x());
    EXPECT_EQ(3, v1.y());
}

TEST(Vec2DTest, Norm)
{
    Vec2D v1{3, 4};
    double norm = Vec2D::Norm(v1);
    EXPECT_EQ(5, norm);
}

TEST(Vec2DTest, Normalize)
{
    Vec2D v1{3, 4};
    v1 = Vec2D::Normalize(v1);
    EXPECT_EQ(3.0/5.0, v1.x());
    EXPECT_EQ(4.0/5.0, v1.y());
}


TEST(DomainTest, Constructor)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {false, true, true, true, false};
    std::vector<bool> r2 = {false, false, true, false, false};
    v.push_back(std::vector<bool>(5, false));
    v.push_back(r2);
    v.push_back(r1);
    v.push_back(r2);
    v.push_back(std::vector<bool>(5, false));

    DomainTest testdomain{v};

    EXPECT_EQ(1.0/4.0, testdomain.domain_.GetH());

    //testdomain.domain_.Dump(std::cout);
}

TEST(BoundaryTest, constructor)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {false, true, true, true, false};
    std::vector<bool> r2 = {false, false, true, true, false};
    v.push_back(std::vector<bool>(5, false));
    v.push_back(r2);
    v.push_back(r1);
    v.push_back(r2);
    v.push_back(std::vector<bool>(5, false));

    BoundaryTest test(v, 1.0/4.0);
    //test.boundary_.Dump(std::cout);
}

TEST(BoundaryTest, IsBoundary)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {true, false, false};
    v.push_back(r1);
    std::vector<bool> r2 = {false, true, false};
    v.push_back(r2);
    v.push_back(r1);

    BoundaryTest test(v, 1.0/2.0);

    EXPECT_EQ(true, test.IsBoundary(1, 1, v));
}

TEST(BCDataTest, GenerateSize)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {true, false, false};
    v.push_back(r1);
    std::vector<bool> r2 = {false, true, false};
    v.push_back(r2);
    v.push_back(r1);

    DomainTest testdomain{v};

    BoundaryDataGeneratorTest testgenerator(testfunction);
    std::vector<double> bc = testgenerator.generator_.Generate(testdomain.domain_);

    srand (time(NULL));
    EXPECT_EQ(3, bc.size());
}



TEST(BCDataTest, GenerateValue)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {true, false, false};
    v.push_back(r1);
    std::vector<bool> r2 = {false, true, false};
    v.push_back(r2);
    v.push_back(r1);

    DomainTest testdomain{v};

    BoundaryDataGeneratorTest testgenerator(testfunction);
    std::vector<double> bc = testgenerator.generator_.Generate(testdomain.domain_);

    srand (time(NULL));
    EXPECT_EQ(1, bc[rand() % bc.size()]);
    EXPECT_EQ(3, bc.size());
}



