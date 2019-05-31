#ifndef SOLVEREIGEN_H
#define SOLVEREIGEN_H


#include "Solver.h"
#include <Eigen/Dense>
#include <complex>

typedef std::complex<double> complex;
typedef Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Vector<complex, Eigen::Dynamic> Vector;

///@brief solver class using Eigen as SoE solver.
class SolverEigen : public Solver
{
    friend class SolverEigenTest;

public:
    virtual std::vector<double> Solve(const Domain& dom) const;

protected:
    ///@brief Function that sets up the voefficient Matrix A.
    Matrix SetupA(const Domain& dom) const;

    ///@brief Function that sets up the RHS of the SoE.
    Vector SetupRhs(const Domain& dom) const;

    ///@brief Setup the matrix necessary to evaluate the solution inside the domain.
    Matrix SetupEval(const Domain& dom) const;

    ///@brief Compute the coefficient R as in Kress (2013)
    double ComputeR(double tj, double t, int n) const;

    ///@brief Kernel as in Kress (2013)
    complex K(int i, int j, const Domain& dom) const;

    ///@brief Kernel as in Kress (2013)
    complex K1(int i, int j, const Domain& dom) const;

    ///@brief Kernel as in Kress (2013)
    complex K2(int i, int j, const Domain& dom) const;

    ///@brief Kernel as in Kress (2013)
    complex Lxy(double x, double y, int i, const Domain& dom) const;

protected:
    
};


#endif /* SOLVEREIGEN_H */