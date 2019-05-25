#include "Domain.h"

#include<algorithm>


Domain::Domain(std::vector<std::vector<bool> > points) :
    points_(points), 
    interior_(GenerateInterior()),
    boundary_(GenerateBoundary()),
    normals_(GenerateNormals())
{
}

std::vector<std::pair<int, int> > Domain::GenerateInterior() const
{
    std::vector<std::pair< int, int> > result; 

    for (int i = 0; i < points_.size(); i++)
    {
        for (int j = 0; j < points_[i].size(); j++)
        {
            if (points_[i][j])
                result.push_back(std::make_pair(i, j));
        }
    }
    std::sort(result.begin(), result.end()); //interior is sorted first by x then by y coordinate

    return result;
}

std::vector<int> Domain::GenerateBoundary() const
{
    std::vector<int> boundaries;
    for (int i = 0; i < interior_.size(); i++)
    {
        //ATTENTION: i am now looking only at the neighbours in x and y direction,
        //not the diagonal ones to determine if something is a boundary element
        std::pair<int, int> coordinates = interior_[i];
        bool is_boundary = false;
        if (coordinates.first != 0 && !points_[coordinates.first-1][coordinates.second])
            is_boundary = true;
        
        if (coordinates.first < points_.size() -1 && !points_[coordinates.first+1][coordinates.second])
            is_boundary = true;

        if (coordinates.second != 0 && !points_[coordinates.first][coordinates.second-1])
            is_boundary = true;
        
        if (coordinates.second < points_[coordinates.first].size() -1 && !points_[coordinates.first][coordinates.second+1])
            is_boundary = true;

        if (is_boundary)
            boundaries.push_back(i);
    }

    return boundaries;
}

std::vector<Vec2D> Domain::GenerateNormals() const
{
    std::vector<Vec2D> normals;
    normals.resize(boundary_.size());

    for (int i = 0; i < boundary_.size(); i++)
    {
        std::pair<int, int> coordinates = interior_[boundary_[i]];
        Vec2D n(0, 0);
        
        if (coordinates.first != 0 && !points_[coordinates.first-1][coordinates.second])
            n += Vec2D(-1, 0);
        
        if (coordinates.first < points_.size() -1 && !points_[coordinates.first+1][coordinates.second])
            n += Vec2D(1, 0);

        if (coordinates.second != 0 && !points_[coordinates.first][coordinates.second-1])
            n += Vec2D(0, -1);
        
        if (coordinates.second < points_[coordinates.first].size() -1 && !points_[coordinates.first][coordinates.second+1])
            n += Vec2D(0, 1);

        n = Vec2D::Normalize(n);

        normals[i] = n;
    }

    return normals;

}



std::ostream& Domain::Dump(std::ostream& ofs) const
{
    for (int i = 0; i < points_.size(); i++)
    {
        for (int j = 0; j < points_[i].size(); j++)
        {
            if (points_[i][j])
                ofs <<  "x";
            else
                ofs << "0";

            ofs << " ";
        }
        ofs << std::endl;
    }

    ofs << std::endl;
    ofs << "Interior" << std::endl;

    for (int i = 0; i < interior_.size(); i++)
    {
        ofs << interior_[i].first << ", " << interior_[i].second << std::endl;
    }

    ofs << std::endl;
    ofs << "Boundary and Normals" << std::endl;

    for (int i = 0; i < boundary_.size(); i++)
    {
        ofs << boundary_[i] << "\t";
        normals_[i].Dump(ofs);
    }
    
    return ofs;    
}
