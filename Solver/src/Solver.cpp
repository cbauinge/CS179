#include "Solver.h"
#include <cmath>
#include <complex>


using namespace std::complex_literals;

//c++ has bessel functions of first and second kind in cmath
//double _j0(double x); first kind order 0
//double _j1(double x); first kind order 1
//double _jn(int n,double x); first kind order n
//double _y0(double x); second kind order 0
//double _y1(double x); second kind order 1
//double _yn(int n,double x); second kind order n

void Solver::SetWaveNumber(double k)
{
    k_ = k;
}

void Solver::SetBoundaryCondition(const std::vector<double>& bc)
{
    bc_ = bc;
}


std::complex<double> Solver::Hankel(int alpha, double val) const
{
    return std::cyl_bessel_j(alpha, val) + 1i*std::cyl_neumann(alpha, val);
}