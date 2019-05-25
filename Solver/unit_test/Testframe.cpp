#include "Solver.h"

///NOTE: this construction is necessary to have access to private members of the tested class.


class SolverTest : public ::testing::Test
{
public:
    virtual void SetUp(){}
    virtual void TearDown(){}
    virtual void TestBody(){}

    Solver solver_;

};




TEST{SolverTest, Constructor}
{
    SolverTest test;
}