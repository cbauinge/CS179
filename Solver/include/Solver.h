#ifndef SOLVER_H
#define SOLVER_H

#include "Domain.h"
#include <vector>


class Solver
{
public:

    ///Solve routine which takes a domain and computes the solution of the Hemholtz equation for all interior 
    ///points. Returns them in the same order as in the domain.
    std::vector<double> Solve(const Domain& dom) const;

    void SetWaveNumber(double k) {k_ = k;}
    void SetBoundaryCondition(const std::vector<double>& bc) {bc_ = bc;}


private:
    double k_; //wave number
    std::vector<double> bc_; //boundary condition
};



#endif /* SOLVER_H */