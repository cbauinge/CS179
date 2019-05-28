#ifndef SOLVER_H
#define SOLVER_H

#include "Domain.h"
#include <vector>
#include <complex>


class Solver
{
    friend class SolverTest;

public:

    ///Solve routine which takes a domain and computes the solution of the Hemholtz equation for all interior 
    ///points. Returns them in the same order as in the domain.
    virtual std::vector<double> Solve(const Domain& dom) const = 0;

    void SetWaveNumber(double k);
    void SetBoundaryCondition(const std::vector<double>& bc);    

protected:
    ///Hankel function of first kind
    std::complex<double> Hankel(int alpha, double val) const;


protected:
    double k_; //wave number
    std::vector<double> bc_; //boundary condition
};



#endif /* SOLVER_H */