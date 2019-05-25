#include "BoundaryDataGenerator.h"


BoundaryDataGenerator::BoundaryDataGenerator(double (*func)(double, double)) :
    bcfunc(func)
{

}

std::vector<double> BoundaryDataGenerator::Generate(const Domain& dom)
{
    std::vector<double> bc;
    bc.resize(dom.GetBoundary().size());
    double h = dom.GetH();

    for (int i = 0; i < bc.size(); i++)
    {
        std::pair<int, int> coordinates = dom.GetInterior()[dom.GetBoundary()[i]];
        bc[i] = bcfunc(h*coordinates.first, h*coordinates.second);
    }

    return bc;
}