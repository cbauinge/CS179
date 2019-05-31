#include <cstdlib>
#include <vector>
#include <exception>
#include <string>
#include <iostream>

#include "Domain.h"
#include "BoundaryDataGenerator.h"
#include "Reader.h"
#include "Writer.h"
#include "Solver.h"
#include "SolverEigen.h"

double bcdata(double x, double y)
{
    return x;
}

int main(int argc, char * argv [])
{
    if (argc != 3) {
        std::cout << "Usage: ./CS179 <path-to-inputfile> <wave number>" << std::endl;
        return 1;
    }

    try
    {
        //set the wave number
        double k = std::stod(argv[2]);

        //Get the domain from the file
        Reader r;
        Domain dom = r.Read(argv[1]);

        //Generate bc according to the above function
        BoundaryDataGenerator bc_generator(bcdata);
        std::vector<double> bc = bc_generator.Generate(dom);

        //Solver the equation with the given wave number and boundary condition
        Solver* solver = new SolverEigen;
        solver->SetWaveNumber(k);
        solver->SetBoundaryCondition(bc);
        std::vector<double> result = solver->Solve(dom);
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