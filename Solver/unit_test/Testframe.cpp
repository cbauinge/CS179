#include "Solver.h"
#include "SolverEigen.h"
#include "gtest/gtest.h"

///NOTE: this construction is necessary to have access to private members of the tested class.


class SolverEigenTest : public ::testing::Test
{
public:
    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Matrix SetupR(const Domain& dom) const {return solver_.SetupR(dom);}

    SolverEigen solver_;

};




TEST(SolverEigenTest, Constructor)
{
    SolverEigenTest test;
}


TEST(SolverEigenTest, SetupR)
{
    std::vector<std::vector<bool> > v;
    std::vector<bool> r1 = {false, true, true, true, false};
    std::vector<bool> r2 = {false, false, true, false, false};
    v.push_back(std::vector<bool>(5, false));
    v.push_back(r2);
    v.push_back(r1);
    v.push_back(r2);
    v.push_back(std::vector<bool>(5, false));

    Domain dom(v);

    SolverEigenTest test;

    Matrix R = test.SetupR(dom);

    std::cout << R << std::endl;
}