#include "BoundaryDataGenerator.h"


BoundaryDataGenerator::BoundaryDataGenerator(double (*func)(double, double)) :
    bcfunc(func)
{

}

std::vector<double> BoundaryDataGenerator::Generate(const Domain& dom)
{
    std::vector<double> bc;
    bc.resize(dom.GetBoundary().size());

    for (int i = 0; i < bc.size(); i++)
    {
        bc[i] = bcfunc(dom.GetBoundary()[i].x, dom.GetBoundary()[i].y);
    }

    return bc;
}