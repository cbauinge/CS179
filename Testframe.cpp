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
    #include "SolverGPU.h"
#endif /*USE_CUDA*/

double bcdata(double x, double y)
{
    return 10*x -5;
}

int main(int argc, char * argv [])
{
    if (argc != 1) 
    {
        std::cout << "Usage: ./Test" << std::endl;
        return 1;
    }

    try
    {
        //set the wave number
        int k = 25;
        //inititalize the boundary data generator
        BoundaryDataGenerator bc_generator(bcdata);

        //Get the domain from the file for different cases
        for (int i = 1; i < 6; i++)
        {
            std::string filename(std::string("../Testdata/test") + std::string(std::to_string(i)) + std::string(".txt"));
            std::cout << "Start reading file " << filename << "..." << std::endl;
            Reader r;
            Domain dom = r.Read(filename.c_str());
            std::cout << "...Finished reading file." << std::endl;
            std::cout << "Domain consists of " << dom.GetBoundary().size() << " boundary elements and "
                << dom.GetInteriorWOBoundary().size() << " interior elements." << std::endl;

            //Generate bc according to the above function
            std::vector<double> bc = bc_generator.Generate(dom);

            //Solver the equation with the given wave number and boundary condition
            #ifdef USE_CUDA
            {
                std::cout << "Start GPU solving process..." << std::endl;
                Solver* solver = new SolverGPU;
                solver->SetWaveNumber(k);
                solver->SetBoundaryCondition(bc);
                auto start = std::chrono::high_resolution_clock::now();
                std::vector<double> result = solver->Solve(dom);
                auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
                std::cout << "...Finished GPU solver. Solver took " << diff.count()/1000.0 << " seconds." << std::endl;
                delete solver;

                //Write the result to the harddrive
                Writer writer;
                std::string outfile = std::string("result_GPU_") + std::string(std::to_string(i)) + std::string(".csv");
                writer.Write(outfile.c_str(), dom, result, bc);
                std::cout << std::endl;
            }
            #endif

            std::cout << "Start CPU solving process..." << std::endl;
            Solver* solver = new SolverEigen;          
            solver->SetWaveNumber(k);
            solver->SetBoundaryCondition(bc);
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<double> result = solver->Solve(dom);
            auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
            std::cout << "...Finished CPU solver. Solver took " << diff.count()/1000.0 << " seconds." << std::endl;
            delete solver;

            //Write the result to the harddrive
            Writer writer;
            std::string outfile = std::string("result_CPU_") + std::string(std::to_string(i)) + std::string(".csv");
            writer.Write(outfile.c_str(), dom, result, bc);

            std::cout << std::endl << std::endl;
        }
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 2;
    }


    return 0;
}