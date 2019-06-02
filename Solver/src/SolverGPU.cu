#include "SolverGPU.cuh"

#define PI 3.14159265359

///@brief Compute the coefficient R as in Kress (2013)
__device__ double ComputeR(double tj, double t, int fn)
{
    double sum = 0.0;

    for (int m = 1; m < fn; m++)
    {
        sum += 1.0/m * cos(m*(t - tj)) + 1.0/(2.0*fn)*cos(fn*(t - tj)); 
    }
    sum /= -fn;

    return sum;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi)
{
    cuDoubleComplex val;

    if (abs(xi-xj) + abs(yi - yj) < 1e-10)
    {
        val = make_cuDoubleComplex(-1.0/(2.0*PI*(dxi*dxi + dyi*dyi)) * (dxi*ddyi - dyi*ddxi), 0.0);
    }
    else
    {
        double diffxji = xj - xi;
        double diffyji = yj - yi;
        double normxji = sqrt(diffxji*diffxji + diffyji*diffyji);
        val = cuCmul(make_cuDoubleComplex(0.0, -kglob/(2.0*normxji)*(dyj*diffxji - dxj*diffyji)), hankel11(kglob*normxji));
    }

    return val;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K1(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi)
{
    cuDoubleComplex val = make_cuDoubleComplex(0.0, 0.0);

    if (abs(xi-xj) + abs(yi - yj) > 1e-10)
    {
        double diffxij = xi - xj;
        double diffyij = yi - yj;
        double normxij = sqrt(diffxij*diffxij + diffyij*diffyij);
        val =  make_cuDoubleComplex(kglob/(2.0*PI*normxij)*(dyj*diffxij - dxj*diffyij)*j1(kglob*normxij), 0.0);
    }

    return val;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K2(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double dti, double dtj)
{
    cuDoubleComplex val = make_cuDoubleComplex(0.0, 0.0);

    if (abs(xi-xj) + abs(yi - yj) > 1e-10)
    {
        double tmp = sin((dti - dtj)/2.0);
        val = cuCsub(K(xi, yi, xj, yj, dxi, dyi, dxj, dyj, ddxi, ddyi),
            cuCmul(K1(xi, yi, xj, yj, dxi, dyi, dxj, dyj, ddxi, ddyi), make_cuDoubleComplex(log(4.0*tmp*tmp), 0.0)));
    }
    else
    {
        val = K(xi, yi, xj, yj, dxi, dyi, dxj, dyj, ddxi, ddyi);
    }

    return val;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex Lxy(double x, double y, double xi, double yi, double dxi, double dyi)
{
    double distx = x - xi;
    double disty = y - yi;
    double normdist = sqrt(distx*distx + disty*disty);
    return cuCmul(make_cuDoubleComplex(0.0, kglob * dtglob /(4.0*normdist)*(dyi*distx - dxi*disty)), hankel11(kglob*normdist));
}

__device__ cuDoubleComplex hankel11(double x)
{
    return make_cuDoubleComplex(j1(x), y1(x));
}


void SetupAKernelCall(int nr_blocks, int nr_threads, cuDoubleComplex* A, 
    double* posx,
    double* posy,
    double* dx,
    double* dy,
    double* ddx,
    double* ddy)
{
    SetupAKernel<<<nr_blocks, nr_threads>>>(A, posx, posy, dx, dy, ddx, ddy);
}

__global__ void SetupAKernel(cuDoubleComplex* A, 
    double* posx,
    double* posy,
    double* dx,
    double* dy,
    double* ddx,
    double* ddy)
{
    uint thread_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (thread_index >= nglob*nglob)
        return;

    int i = thread_index / nglob;
    int j = thread_index % nglob;

    cuDoubleComplex r = make_cuDoubleComplex(ComputeR(dtglob*j, dtglob*i, nglob/2), 0.0);
    cuDoubleComplex k1 = K1(posx[i], posy[i], posx[j], posy[j], dx[i], dy[i], dx[j], dy[j], ddx[i], ddy[i]);
    cuDoubleComplex k2 = K2(posx[i], posy[i], posx[j], posy[j], dx[i], dy[i], dx[j], dy[j], ddx[i], ddy[i], dtglob*i, dtglob*j);

    cuDoubleComplex a = cuCadd(cuCmul(r, k1), cuCmul(make_cuDoubleComplex(dtglob, 0.0),k2));
    A[thread_index] = a;
    //int fn = floor(n/2);

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         A(i, j) = ComputeR(dt*j, dt*i, fn)*K1(i, j, dom) + dt*K2(i, j, dom);
    //     }
    // }
}

__global__ void ComputeIdmAKernel(cuDoubleComplex* A)
{
    uint thread_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (thread_index >= nglob*nglob)
        return;

    int i = thread_index / nglob;
    int j = thread_index % nglob;

    cuDoubleComplex a = A[thread_index];
    cuDoubleComplex id = i == j ? make_cuDoubleComplex(1.0, 0.0) :
        make_cuDoubleComplex(0.0, 0.0);
    
    A[thread_index] = cuCsub(id, a);
}

void ComputeIdmAKernelCall(int nr_blocks, int nr_threads, cuDoubleComplex* A)
{
    ComputeIdmAKernel<<<nr_blocks, nr_threads>>>(A);
}


__global__ void SetupRhsKernel(cuDoubleComplex* rhs)
{

}

void SetupEvalKernelCall(int nr_blocks, int nr_threads, cuDoubleComplex* Eval,
    double* posx_int,
    double* posy_int,
    double* posx_bound,
    double* posy_bound,
    double* dx,
    double* dy,
    int m)
{
    SetupEvalKernel<<<nr_blocks, nr_threads>>>(Eval, 
        posx_int, 
        posy_int, 
        posx_bound, 
        posy_bound,
        dx,
        dy,
        m);
}


__global__ void SetupEvalKernel(cuDoubleComplex* Eval,
    double* posx_int,
    double* posy_int,
    double* posx_bound,
    double* posy_bound,
    double* dx,
    double* dy,
    int m)
{

    uint thread_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (thread_index >= m*nglob)
        return;

    // int i = thread_index / nglob;
    // int j = thread_index % nglob;

    int i = thread_index % m;
    int j = thread_index / m;

    Eval[thread_index] = Lxy(posx_int[i], posy_int[i], posx_bound[j], posy_bound[j], dx[j], dy[j]);
}


void copyToConstMemory(double k_, double dt_, int n_)
{
    checkCudaErrors(cudaMemcpyToSymbol(kglob, &k_, sizeof(double)));
    checkCudaErrors(cudaMemcpyToSymbol(dtglob, &dt_, sizeof(double)));
    checkCudaErrors(cudaMemcpyToSymbol(nglob, &n_, sizeof(int)));

}