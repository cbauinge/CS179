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
    size_t freemem, totalmem;
    checkCudaErrors(cudaMemGetInfo	(	&freemem, &totalmem));
    std::cout << "Free mem (mb) = " << freemem/1024/1024 << " Total mem (mb)= " << totalmem/1024/1024 << std::endl;

    std::cout << "Minimum necesary memory (mb)= " << m*n*sizeof(cuDoubleComplex)/1024/1024 << std::endl;

    double* devx, *devy, *devdx, *devdy, *devddx, *devddy; //data on dpu
    SetupGPUData(dom, &devx, &devy, &devdx, &devdy, &devddx, &devddy);

    cuDoubleComplex* A;
    SetupA(&A, devx, devy, devdx, devdy, devddx, devddy, dom.GetBoundary().size());
    {
        std::ofstream ofs;
        cuDoubleComplex* host_A = (cuDoubleComplex*) malloc(n*n*sizeof(cuDoubleComplex));
        checkCudaErrors(cudaMemcpy(host_A, A, n*n*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
        ofs.open("TestMatrix.txt");
        Dump(ofs, host_A, n, n);
        ofs.close();
        free(host_A);
    }
    //I-A stored in A
    ComputeIdmA(A, n);
    {
        std::ofstream ofs;
        cuDoubleComplex* host_A = (cuDoubleComplex*) malloc(n*n*sizeof(cuDoubleComplex));
        checkCudaErrors(cudaMemcpy(host_A, A, n*n*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
        ofs.open("TestIdmAMatrix.txt");
        Dump(ofs, host_A, n, n);
        ofs.close();
        free(host_A);
    }
    
    cuDoubleComplex* b;
    SetupRhs(dom, &b);
    {
        std::ofstream ofs;
        cuDoubleComplex* host_b = (cuDoubleComplex*) malloc(n*sizeof(cuDoubleComplex));
        checkCudaErrors(cudaMemcpy(host_b, b, n*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
        ofs.open("TestRhs.txt");
        Dump(ofs, host_b, n);
        ofs.close();
        free(host_b);
    }


    //call cusolver
    int bufferSize = 0;
    int *info = nullptr;
    cuDoubleComplex *buffer = nullptr;
    int *ipiv = NULL;
    int lda = n;
    cusolverDnHandle_t handle;
    checkCudaErrors(cusolverDnCreate(&handle));
    checkCudaErrors(cusolverDnZgetrf_bufferSize(handle, n, n, A, lda, &bufferSize));
    checkCudaErrors(cudaMalloc(&info, sizeof(int)));
    checkCudaErrors(cudaMalloc(&buffer, sizeof(cuDoubleComplex)*bufferSize));
    //checkCudaErrors(cudaMalloc(&ipiv, sizeof(int)*n));

    //do the LU factorization
    checkCudaErrors(cusolverDnZgetrf(handle, n, n, A, lda, buffer, ipiv, info));
    {
        int* host_info = (int*) malloc(sizeof(int));
        checkCudaErrors(cudaMemcpy(host_info, info, sizeof(int), cudaMemcpyDeviceToHost));
        if (*host_info != 0)
            std::cout << "something we wrong when LU the SoE on the device" << std::endl;
        free(host_info);
    }

    checkCudaErrors(cusolverDnZgetrs(handle, CUBLAS_OP_T, n, 1, A, lda, ipiv, b, n, info));
    checkCudaErrors(cudaDeviceSynchronize());
    {
        int* host_info = (int*) malloc(sizeof(int));
        checkCudaErrors(cudaMemcpy(host_info, info, sizeof(int), cudaMemcpyDeviceToHost));
        if (*host_info != 0)
            std::cout << "something we wrong when solving the SoE on the device" << std::endl;
        free(host_info);
    }

    checkCudaErrors(cusolverDnDestroy(handle));
    checkCudaErrors(cudaFree(ipiv));
    checkCudaErrors(cudaFree(buffer));
    checkCudaErrors(cudaFree(info));
    checkCudaErrors(cudaFree(A));

    {
        cuDoubleComplex* host_b = (cuDoubleComplex*) malloc(n*sizeof(cuDoubleComplex));
        checkCudaErrors(cudaMemcpy(host_b, b, n*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
        std::ofstream ofs;
        ofs.open("DensityTest.txt");
        Dump(ofs, host_b, n);
        ofs.close();
        std::cout << "Dumped density device" << std::endl;
        free(host_b);
    }

    //evaluate result on every interior point
    cuDoubleComplex* Eval;
    SetupEval(dom, &Eval, devx, devy, devdx, devdy);
    checkCudaErrors(cudaDeviceSynchronize());

    cuDoubleComplex* hostEval = (cuDoubleComplex*) malloc(m*n*sizeof(cuDoubleComplex));
    checkCudaErrors(cudaMemcpy(hostEval, Eval, m*n*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    {
        std::ofstream ofs;
        ofs.open("TestEvalMatrix.txt");
        Dump(ofs, hostEval, m, n);
        ofs.close();
    }

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

    checkCudaErrors(cublasDestroy(cublashandle));
    checkCudaErrors(cudaFree(Eval));
    checkCudaErrors(cudaFree(b));
    checkCudaErrors(cudaDeviceSynchronize());

    //copy the result to result vector.
    cuDoubleComplex* output = (cuDoubleComplex*) malloc(m*sizeof(cuDoubleComplex));
    checkCudaErrors(cudaMemcpy(output,doutput, m*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(doutput));

    std::vector<double> result(m);
    #pragma omp parallel for
    for (int k = 0; k < m; k++)
    {
        result[k] = cuCreal(output[k]);
    }

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