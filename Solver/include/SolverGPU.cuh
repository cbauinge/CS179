#include "SolverGPU.h"

__device__ __constant__ double kglob;
__device__ __constant__ double dtglob;
__device__ __constant__ int nglob;

///@brief Compute the coefficient R as in Kress (2013)
__device__ double ComputeR(double tj, double t, int fn);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K1(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K2(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double dti, double dtj);

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex Lxy(double x, double y, double xi, double yi, double dxi, double dyi);

__device__ cuDoubleComplex hankel11(double x);

__global__ void SetupAKernel(cuDoubleComplex* A, 
    double* posx,
    double* posy,
    double* dx,
    double* dy,
    double* ddx,
    double* ddy);

void SetupAKernelCall(int nr_blocks, int nr_threads, cuDoubleComplex* A, 
    double* posx,
    double* posy,
    double* dx,
    double* dy,
    double* ddx,
    double* ddy);

__global__ void ComputeIdmAKernel(cuDoubleComplex* A);

void ComputeIdmAKernelCall(int nr_blocks, int nr_threads, cuDoubleComplex* A);

__global__ void SetupRhsKernel(cuDoubleComplex* rhs);

void SetupEvalKernelCall(int nr_blocks, int nr_threads, cuDoubleComplex* Eval,
    double* posx_int,
    double* posy_int,
    double* posx_bound,
    double* posy_bound,
    double* dx,
    double* dy,
    int m);

__global__ void SetupEvalKernel(cuDoubleComplex* Eval,
    double* posx_int,
    double* posy_int,
    double* posx_bound,
    double* posy_bound,
    double* dx,
    double* dy,
    int m);

void copyToConstMemory(double k_, double dt_, int n_);