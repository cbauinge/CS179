#include "SolverGPU.h"

SolverGPU::SolverGPU()
{

}


std::vector<double> SolverGPU::Solve(const Domain& dom) const
{
    std::vector<double> result;

    cuDoubleComplex* A;
    SetupA(dom, A);

    cuDoubleComplex* b;
    SetupRhs(dom, b);

    //call cusolver

    //evaluate result on every interior point
    cuDoubleComplex* Eval;
    SetupEval(dom, Eval);

    //call cusolver to multiply eval*the result from the cusolver call;

    //copy the result to result vector.


    checkCudaErrors(cudaFree(A));
    checkCudaErrors(cudaFree(b));
    checkCudaErrors(cudaFree(Eval));

    return result;
}


void SolverGPU::SetupA(const Domain& dom, cuDoubleComplex* A) const
{

}

void SolverGPU::SetupRhs(const Domain& dom, cuDoubleComplex* rhs) const
{
    // int n = dom.GetBoundary().size();
    // Vector b(n);

    // for (int i = 0; i < n; i++)
    // {
    //     b(i) = -2.0*bc_[i];
    // }
}

void SolverGPU::SetupEval(const Domain& dom, cuDoubleComplex* eval) const
{
 //Matrix Eval(dom.GetInteriorWOBoundary().size(), dom.GetBoundary().size());
}