#include "SolverGPU.h"


///@brief Compute the coefficient R as in Kress (2013)
__device__ double ComputeR(double tj, double t, int n);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K(int i, int j, const Domain& dom);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K1(int i, int j, const Domain& dom);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K2(int i, int j, const Domain& dom);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex Lxy(double x, double y, int i, const Domain& dom);

__device__ cuDoubleComplex hankel11(double x);


__global__ void SetupAKernel(cuDoubleComplex* A);

__global__ void SetupRhsKernel(cuDoubleComplex* rhs);

__global__ void SetupEvalKernel(cuDoubleComplex* Eval);