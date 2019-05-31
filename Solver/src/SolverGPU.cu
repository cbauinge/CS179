#include "SolverGPU.cuh"

#define PI 3.14159265359

///@brief Compute the coefficient R as in Kress (2013)
__device__ double ComputeR(double tj, double t, int n)
{
    double sum = 0.0;

    for (int m = 1; m < n; m++)
    {
        sum += 1.0/m * cos(m*(t - tj)) + 1/(2.0*n)*cos(n*(t - tj)); 
    }
    sum /= -n;

    return sum;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double k)
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
        val = cuCmul(make_cuDoubleComplex(0.0, -k/(2.0*normxji)*(dyj*diffxji - dxj*diffyji)), hankel11(k*normxji));
    }

    return val;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K1(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double k)
{
    cuDoubleComplex val = make_cuDoubleComplex(0.0, 0.0);

    if (abs(xi-xj) + abs(yi - yj) > 1e-10)
    {
        double diffxij = xi - xj;
        double diffyij = yi - yj;
        double normxij = sqrt(diffxij*diffxij + diffyij*diffyij);
        val =  make_cuDoubleComplex(k/(2.0*PI*normxij)*(dyj*diffxij - dxj*diffyij)*j1(k*normxij), 0.0);
    }

    return val;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex K2(double xi, double yi, double xj, double yj, double dxi, 
    double dyi, double dxj, double dyj, double ddxi, double ddyi, double dti, double dtj, double k)
{
    cuDoubleComplex val = make_cuDoubleComplex(0.0, 0.0);

    if (abs(xi-xj) + abs(yi - yj) > 1e-10)
    {
        double tmp = sin((dti - dtj)/2.0);
        val = cuCsub(K(xi, yi, xj, yj, dxi, dyi, dxj, dyj, ddxi, ddyi, k),
            cuCmul(K1(xi, yi, xj, yj, dxi, dyi, dxj, dyj, ddxi, ddyi, k), make_cuDoubleComplex(log(4.0*tmp*tmp), 0.0)));
    }
    else
    {
        val = K(xi, yi, xj, yj, dxi, dyi, dxj, dyj, ddxi, ddyi, k);
    }

    return val;
}

///@brief Kernel as in Kress (2013)
__device__ cuDoubleComplex Lxy(double x, double y, double xi, double yi, double dxi, double dyi, double k, double dt)
{
    double distx = x - xi;
    double disty = y - yi;
    double normdist = sqrt(distx*distx + disty*disty);
    return cuCmul(make_cuDoubleComplex(0.0, k * dt /(4.0*normdist)*(dyi*distx - dxi*disty)), hankel11(k*normdist));
}

__device__ cuDoubleComplex hankel11(double x)
{
    return make_cuDoubleComplex(j1(x), y1(x));
}


__global__ void SetupAKernel(cuDoubleComplex* A, int n, double dt)
{
    //int fn = floor(n/2);

    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         A(i, j) = ComputeR(dt*j, dt*i, fn)*K1(i, j, dom) + dt*K2(i, j, dom);
    //     }
    // }
}

__global__ void SetupRhsKernel(cuDoubleComplex* rhs)
{

}


__global__ void SetupEvalKernel(cuDoubleComplex* Eval)
{
    // for (int i = 0; i < Eval.rows(); i++)
    // {
    //     for (int j = 0; j < Eval.cols(); j++)
    //     {   
    //         std::pair<int, int> coordinates = dom.GetInteriorWOBoundary()[i];
    //         Eval(i, j) = Lxy(coordinates.first*dom.GetH(), coordinates.second*dom.GetH(), j, dom);
    //     }
    // }
}