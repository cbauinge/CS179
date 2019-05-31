#ifndef SOLVERGPU_H
#define SOLVERGPU_H

#ifdef USE_CUDA


#include "Solver.h"
#include <complex>
#include <vector>
#include "cusolverDn.h"
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "helper_cuda.h"


#define checkCudaErr


///@brief solver class using cusolver as SoE solver.
class SolverGPU : public Solver
{
    friend class SolverGPUTest;

public:
    SolverGPU();

    virtual std::vector<double> Solve(const Domain& dom) const;

protected:
    ///@brief Function that sets up the voefficient Matrix A.
    ///@param [inout] A. Assumes it is not allocated, allocates it.
    void SetupA(const Domain& dom, cuDoubleComplex* A) const;

    ///@brief Function that sets up the RHS of the SoE.
    ///@param [inout] rhs: assumes rhs to not be allocated. allocates it in the function.
    void SetupRhs(const Domain& dom, cuDoubleComplex* rhs) const;

    ///@brief Setup the matrix necessary to evaluate the solution inside the domain.
    void SetupEval(const Domain& dom, cuDoubleComplex* eval) const;
protected:
    
};

#endif /* USE_CUDA */


#endif /* SOLVERGPU_H */