#include "SolverGPU.h"


///@brief Compute the coefficient R as in Kress (2013)
__device__ double ComputeR(double tj, double t, int n);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double k);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K1(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double k);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K2(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double dti, double dtj, double k);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex Lxy(double x, double y, double xi, double yi, double dxi, double dyi, double k, double dt);

__device__ cuDoubleComplex hankel11(double x);


__global__ void SetupAKernel(cuDoubleComplex* A, 
    double* posx,
    double* posy,
    double* dx,
    double* dy,
    double* ddx,
    double* ddy,
    int n, double dt, double k);

__global__ void SetupRhsKernel(cuDoubleComplex* rhs);

__global__ void SetupEvalKernel(cuDoubleComplex* Eval,
    double* posx_int,
    double* posy_int,
    double* posx_bound,
    double* posy_bound,
    double* dx,
    double* dy,
    int m, 
    int n,
    double k,
    double dt);