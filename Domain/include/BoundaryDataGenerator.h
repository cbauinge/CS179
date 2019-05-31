#ifndef BOUNDARYDATAGENERATOR_H
#define BOUNDARYDATAGENERATOR_H

#include "Domain.h"
#include <vector>


/// @brief Class that generates the boundary condition for each point on the boundary
/// according to the given function.
class BoundaryDataGenerator
{
    friend class BoundaryDataGeneratorTest;

public:
    BoundaryDataGenerator(double (*func)(double, double));

    std::vector<double> Generate(const Domain& dom);
private:
    double (*bcfunc)(double, double);
};

#endif /* BOUNDARYDATAGENERATOR_H */

