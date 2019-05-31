#ifndef SOLVER_H
#define SOLVER_H

#include "Domain.h"
#include <vector>
#include <complex>

/// @brief Abstract solver class. 
/// Every solver should be derived from this class. There will be own child class
/// for the CPU version using EIGEN and one child class for the GPU version using
/// cusolver.
class Solver
{
    friend class SolverTest;

public:

    ///@brief Solve routine which takes a domain and computes the solution of the Hemholtz equation for all interior 
    ///points. 
    ///@return solution at interior points of the domain in same order as in the domain.
    virtual std::vector<double> Solve(const Domain& dom) const = 0;

    ///@brief Setter for the wave_number.
    void SetWaveNumber(double k);
    
    ///@brief Setter for the boundary condition.
    void SetBoundaryCondition(const std::vector<double>& bc);    

protected:
    ///@brief Hankel function of first kind.
    std::complex<double> Hankel(int alpha, double val) const;


protected:
    double k_; //</wave number.
    std::vector<double> bc_; ///<boundary condition ordered as the boundary condition in the domain.
};



#endif /* SOLVER_H */