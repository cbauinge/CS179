#include "Solver.h"
#include <cmath>
#include <complex>


//c++ has bessel functions of first and second kind in cmath
//double _j0(double x); first kind order 0
//double _j1(double x); first kind order 1
//double _jn(int n,double x); first kind order n
//double _y0(double x); second kind order 0
//double _y1(double x); second kind order 1
//double _yn(int n,double x); second kind order n

std::vector<double> Solver::Solve(const Domain& dom) const
{
    std::vector<double> result;
    result.resize(dom.GetInterior().size()); //each interior point gets a result.

    //TODO: setup the coefficient matrix -> need function that integrates over all the boundary values.
    //TODO: solve the coefficient matrix
    //TODO: compute the solution for all the interior points and save them in result.
    
    return result;
}