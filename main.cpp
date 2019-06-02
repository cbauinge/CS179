#include <cstdlib>
#include <vector>
#include <exception>
#include <string>
#include <iostream>
#include <chrono>

#include "Domain.h"
#include "BoundaryDataGenerator.h"
#include "Reader.h"
#include "Writer.h"
#include "Solver.h"
#include "SolverEigen.h"

#ifdef USE_CUDA
    //#undef USE_CUDA
    #include "SolverGPU.h"
#endif /*USE_CUDA*/

double bcdata(double x, double y)
{
    return 10*x -5;
}

int main(int argc, char * argv [])
{
    if (argc != 3 && argc != 1) 
    {
        std::cout << "Usage: ./CS179 for default OR ./CS179 <path-to-inputfile> <wave number>" << std::endl;
        return 1;
    }

    try
    {
        //set the wave number
        double k = argc == 3 ? std::stod(argv[2]) : 25;

        //Get the domain from the file
        Reader r;
        Domain dom = argc == 3 ? r.Read(argv[1]) : r.Read("../Testdata/test.txt");

        //Generate bc according to the above function
        BoundaryDataGenerator bc_generator(bcdata);
        std::vector<double> bc = bc_generator.Generate(dom);

        //Solver the equation with the given wave number and boundary condition
        #ifdef USE_CUDA
            Solver* solver = new SolverGPU;
        #else
            Solver* solver = new SolverEigen;
        #endif
        
        solver->SetWaveNumber(k);
        solver->SetBoundaryCondition(bc);
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> result = solver->Solve(dom);
        auto diff = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start);
        std::cout << "Solver took " << diff.count() << " seconds." << std::endl;
        delete solver;

        //Write the result to the harddrive
        Writer writer;
        writer.Write("result.csv", dom, result, bc);
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 2;
    }


    return 0;
}