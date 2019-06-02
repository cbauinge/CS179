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
#include <fstream>
#include <iomanip>


///@brief solver class using cusolver as SoE solver.
class SolverGPU : public Solver
{
    friend class SolverGPUTest;

public:
    SolverGPU();

    virtual std::vector<double> Solve(const Domain& dom) const;

protected:
    ///@brief Function that sets up the voefficient Matrix A!
    ///@param [inout] A. Assumes it is not allocated, allocates it.
    void SetupA(cuDoubleComplex** A, 
        double *devx, 
        double *devy,
        double *devdx,
        double *devdy,
        double *devddx, 
        double *devddy,
        int n) const;

    ///@brief Compute I-A
    void ComputeIdmA(cuDoubleComplex* A, int n) const;

    ///@brief Function that sets up the RHS of the SoE.
    ///@param [inout] rhs: assumes rhs to not be allocated. allocates it in the function.
    void SetupRhs(const Domain& dom, cuDoubleComplex** rhs) const;

    ///@brief Setup the matrix necessary to evaluate the solution inside the domain.
    void SetupEval(const Domain& dom, cuDoubleComplex** Eval,
        double* posx_bound,
        double* posy_bound,
        double* dx,
        double* dy) const;

    ///@brief Get the data from the boundary in array form and move them to the GPU
    void SetupGPUData(const Domain& dom, 
        double **devx, 
        double **devy, 
        double **devdx, 
        double **devdy,
        double **devddx, 
        double **devddy) const;
protected:
    ///@brief small helper function that writes the real part of a matrix
    ///A to the given ostream.
    void Dump(std::ofstream& ofs, cuDoubleComplex* A, int m, int n) const;
    
    ///@brief small helper function that writes a vector to the given ofstream
    void Dump(std::ofstream& ofs, cuDoubleComplex* v, int n) const;

    ///@brief helper function tht checks the cusolver info and outputs problems
    ///@param [in] dinfo device info
    void CheckInfo(int *dinfo) const;

};

#endif /* USE_CUDA */


#endif /* SOLVERGPU_H */