#include "Domain.h"


Domain::Domain(std::vector<std::vector<bool> > points) :
    points_(points), 
    interior_(GenerateInterior()),
    boundary_(GenerateBoundary()),
    normals_(GenerateNormals())
{
}

std::vector<std::pair<int, int> > Domain::GenerateInterior() const
{

}

std::vector<int> Domain::GenerateBoundary() const
{

}

std::vector<Vec2D> Domain::GenerateNormals() const
{

}