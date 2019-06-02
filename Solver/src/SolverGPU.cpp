#include "SolverGPU.h"
#include "SolverGPU.cuh"
#include "cublas_v2.h"
#include <iostream>
#include <fstream>

SolverGPU::SolverGPU()
{

}


std::vector<double> SolverGPU::Solve(const Domain& dom) const
{
    int m = dom.GetInteriorWOBoundary().size();
    int n = dom.GetBoundary().size();

    //check the memory available and estimate the minimal memory requirement.
    size_t freemem, totalmem;
    checkCudaErrors(cudaMemGetInfo	(	&freemem, &totalmem));
    std::cout << "Free memmory (mb) = " << freemem/1024/1024 << std::endl;
    std::cout << "Total mem (mb) = " << totalmem/1024/1024 << std::endl;
    std::cout << "Minimum necesary memory (mb)= " << m*n*sizeof(cuDoubleComplex)/1024/1024 << std::endl;

    //setup the boundary data on the gpu
    double* devx, *devy, *devdx, *devdy, *devddx, *devddy; // boundary data on dpu
    SetupGPUData(dom, &devx, &devy, &devdx, &devdy, &devddx, &devddy);

    //Setup A
    cuDoubleComplex* A;
    SetupA(&A, devx, devy, devdx, devdy, devddx, devddy, dom.GetBoundary().size());

    //I-A stored in A
    ComputeIdmA(A, n);
    
    cuDoubleComplex* b;
    SetupRhs(dom, &b);

    //call cusolver
    int bufferSize = 0;
    int *info = nullptr;
    cuDoubleComplex *buffer = nullptr;
    int *ipiv = NULL;
    int lda = n;
    cusolverDnHandle_t handle;
    checkCudaErrors(cusolverDnCreate(&handle));
    checkCudaErrors(cusolverDnZgetrf_bufferSize(handle, n, n, A, lda, &bufferSize)); //Get buffer size
    checkCudaErrors(cudaMalloc(&info, sizeof(int)));
    checkCudaErrors(cudaMalloc(&buffer, sizeof(cuDoubleComplex)*bufferSize));
    //checkCudaErrors(cudaMalloc(&ipiv, sizeof(int)*n));

    //do the LU factorization
    checkCudaErrors(cusolverDnZgetrf(handle, n, n, A, lda, buffer, ipiv, info)); //LU factorization
    checkCudaErrors(cudaDeviceSynchronize());
    CheckInfo(info);

    checkCudaErrors(cusolverDnZgetrs(handle, CUBLAS_OP_T, n, 1, A, lda, ipiv, b, n, info)); //T because the format is a row major but col major is necessary
    checkCudaErrors(cudaDeviceSynchronize());
    CheckInfo(info);

    //destroy and free stuff necessary for solver
    checkCudaErrors(cusolverDnDestroy(handle));
    checkCudaErrors(cudaFree(ipiv));
    checkCudaErrors(cudaFree(buffer));
    checkCudaErrors(cudaFree(info));
    checkCudaErrors(cudaFree(A));

    //evaluate result on every interior point
    cuDoubleComplex* Eval;
    SetupEval(dom, &Eval, devx, devy, devdx, devdy);
    checkCudaErrors(cudaDeviceSynchronize());

    //call cusolver to multiply eval*the result from the cusolver call;
    cublasHandle_t cublashandle;
    checkCudaErrors(cublasCreate(&cublashandle));
    cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex zero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex* doutput;
    checkCudaErrors(cudaMalloc((void**)&doutput, m*sizeof(cuDoubleComplex)));
    checkCudaErrors(cublasZgemv(cublashandle, 
        CUBLAS_OP_N,
        m, n,
        &one,
        Eval, m,
        b, 1,
        &zero,
        doutput, 1));

    //destroy and free stuff necessary for evaluation
    checkCudaErrors(cublasDestroy(cublashandle));
    checkCudaErrors(cudaFree(Eval));
    checkCudaErrors(cudaFree(b));
    checkCudaErrors(cudaDeviceSynchronize());

    //copy the result to result vector.
    cuDoubleComplex* output = (cuDoubleComplex*) malloc(m*sizeof(cuDoubleComplex));
    checkCudaErrors(cudaMemcpy(output,doutput, m*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(doutput));

    //NOTE: this could be performed on the GPU
    std::vector<double> result(m);
    #pragma omp parallel for
    for (int k = 0; k < m; k++)
    {
        result[k] = cuCreal(output[k]);
    }

    //free memory
    free(output);
    checkCudaErrors(cudaFree(devx));
    checkCudaErrors(cudaFree(devy));
    checkCudaErrors(cudaFree(devdx));
    checkCudaErrors(cudaFree(devdy));
    checkCudaErrors(cudaFree(devddx));
    checkCudaErrors(cudaFree(devddy));

    return result;
}


void SolverGPU::SetupA(cuDoubleComplex** A, 
    double *devx, 
    double *devy,
    double *devdx,
    double *devdy,
    double *devddx, 
    double *devddy,
    int n) const
{
    checkCudaErrors(cudaMalloc((void **)A, n*n*sizeof(cuDoubleComplex)));

    int nr_threads = 128;
    int nr_blocks = n*n/nr_threads +1;

    SetupAKernelCall(nr_blocks, nr_threads, *A, devx, devy, devdx, devdy, devddx, devddy);
    checkCudaErrors(cudaDeviceSynchronize());
}

void SolverGPU::ComputeIdmA(cuDoubleComplex* A, int n) const
{
    int nr_threads = 128;
    int nr_blocks = n*n/nr_threads + 1;

    ComputeIdmAKernelCall(nr_blocks, nr_threads, A);
    checkCudaErrors(cudaDeviceSynchronize());
}

void SolverGPU::SetupRhs(const Domain& dom, cuDoubleComplex** rhs) const
{
    int n = dom.GetBoundary().size();
    
    cuDoubleComplex* tmp = (cuDoubleComplex*) malloc(n*sizeof(cuDoubleComplex));
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        tmp[i] = make_cuDoubleComplex(-2.0*bc_[i], 0.0); //i am not sure if I can directly copy double boundary values on the host to cuDoubleComplex values on the device. 
    }

    checkCudaErrors(cudaMalloc((void **)rhs, n*sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaMemcpy(*rhs, tmp, n*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

    free(tmp);
}

void SolverGPU::SetupEval(const Domain& dom, cuDoubleComplex** Eval,
    double* posx_bound,
    double* posy_bound,
    double* dx,
    double* dy) const
{
    int m = dom.GetInteriorWOBoundary().size();
    int n = dom.GetBoundary().size();
    double* x;
    double* y;
    double* devx;
    double* devy;

    dom.GetArrayDataInterior(&x, &y);
    checkCudaErrors(cudaMalloc((void**)&devx, m*sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&devy, m*sizeof(double)));
    checkCudaErrors(cudaMemcpy(devx, x, m*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(devy, y, m*sizeof(double), cudaMemcpyHostToDevice));
    free(x);
    free(y);

    checkCudaErrors(cudaMalloc((void **)Eval, m*n*sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaMemset(*Eval, 0, m*n*sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaDeviceSynchronize());


    int nr_threads = 128;
    int nr_blocks = m*n/nr_threads + 1;
    SetupEvalKernelCall(nr_blocks, nr_threads, *Eval, devx, devy, posx_bound, posy_bound, dx, dy, m);
    checkCudaErrors(cudaDeviceSynchronize());

    checkCudaErrors(cudaFree(devx));
    checkCudaErrors(cudaFree(devy));
    checkCudaErrors(cudaDeviceSynchronize());
}

void SolverGPU::SetupGPUData(const Domain& dom, 
    double **devx, 
    double **devy, 
    double **devdx, 
    double **devdy,
    double **devddx, 
    double **devddy) const
{
    int size = dom.GetBoundary().size();

    checkCudaErrors(cudaMalloc((void **)devx, size*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)devy, size*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)devdx, size*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)devdy, size*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)devddx, size*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)devddy, size*sizeof(double)));

    double* x, *y, *dx, *dy, *ddx, *ddy;
    dom.GetBoundary().GetArrayData(&x, &y, &dx, &dy, &ddx, &ddy);

    checkCudaErrors(cudaMemcpy(*devx, x, size*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(*devy, y, size*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(*devdx, dx, size*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(*devdy, dy, size*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(*devddx, ddx, size*sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(*devddy, ddy, size*sizeof(double), cudaMemcpyHostToDevice));

    double dt_host = 2.0*M_PI/size;
    copyToConstMemory(k_, dt_host, size);

    free(x);
    free(y);
    free(dx);
    free(dy);
    free(ddx);
    free(ddy);
    checkCudaErrors(cudaDeviceSynchronize());
}


void SolverGPU::Dump(std::ofstream& ofs, cuDoubleComplex* A, int m, int n) const
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ofs << "(" << cuCreal(A[j*m + i]) << "," << cuCimag(A[j*m + i]) << ")" << " ";
        }
        ofs << std::endl;
    }
}


void SolverGPU::Dump(std::ofstream& ofs, cuDoubleComplex* v, int n) const
{
    for (int i = 0; i < n; i++)
    {
        ofs << "(" << cuCreal(v[i]) << "," << cuCimag(v[i]) << ")" << std::endl;
    }
}


void SolverGPU::CheckInfo(int *dinfo) const
{
    int* host_info = (int*) malloc(sizeof(int));
    checkCudaErrors(cudaMemcpy(host_info, dinfo, sizeof(int), cudaMemcpyDeviceToHost));
    if (*host_info != 0)
        std::cout << "Sth went wrong, info = " << *host_info << std::endl;
    free(host_info);
}