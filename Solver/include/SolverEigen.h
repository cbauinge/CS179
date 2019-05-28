#ifndef SOLVEREIGEN_H
#define SOLVEREIGEN_H


#include "Solver.h"
#include <Eigen/Dense>
#include <complex>

typedef std::complex<double> complex;
typedef Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Vector<complex, Eigen::Dynamic> Vector;


class SolverEigen : public Solver
{
    friend class SolverEigenTest;

public:
    virtual std::vector<double> Solve(const Domain& dom) const;

protected:
    //Function that sets up the helper Matrix R as in Kress
    Matrix SetupA(const Domain& dom) const;
    Vector SetupRhs(const Domain& dom) const;
    Matrix SetupEval(const Domain& dom) const;

    Matrix SetupR(const Domain& dom) const;
    complex K(int i, int j, const Domain& dom) const;
    complex K1(int i, int j, const Domain& dom) const;
    complex K2(int i, int j, const Domain& dom) const;
    complex Lxy(double x, double y, int i, const Domain& dom) const;

protected:
    
};


#endif /* SOLVEREIGEN_H */