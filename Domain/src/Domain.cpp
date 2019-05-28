#include "Domain.h"

#include<algorithm>


Domain::Domain(std::vector<std::vector<bool> > points) :
    points_(points), 
    interior_(GenerateInterior()),
    interior_without_boundary_(GenerateInteriorWOBoundary()),
    h_(1.0/(points.size()-1.0)), //we scale x direction to 1. from that and the number of points we get h
    boundary_(points, h_)
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


std::vector<std::pair<int, int> > Domain::GenerateInteriorWOBoundary() const
{
    std::vector<std::pair< int, int> > result; 

    for (int i = 0; i < points_.size(); i++)
    {
        for (int j = 0; j < points_[i].size(); j++)
        {
            if (points_[i][j] && !boundary_.IsBoundary(i, j, points_))
                result.push_back(std::make_pair(i, j));
        }
    }
    std::sort(result.begin(), result.end()); //interior is sorted first by x then by y coordinate

    return result;
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
    ofs << "Boundary positions" << std::endl;

    for (int i = 0; i < boundary_.size(); i++)
    {
        ofs << boundary_[i].i << ", " << boundary_[i].j << std::endl;
    }
    
    return ofs;    
}
