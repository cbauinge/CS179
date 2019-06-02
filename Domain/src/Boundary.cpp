#include "Boundary.h"
#include <algorithm>
#include <exception>
#include <cmath>
#include <iostream>


Boundary::Boundary(const std::vector<std::vector<bool> >& points, double h) : 
    nx_(points.size()), ny_(points[0].size())
{
    GeneratePoints(points, h);
    GenerateOrdering(points);
    GenerateGeometricalData(points, h);
}

std::ostream& Boundary::Dump(std::ostream& ofs) const
{
    for (int i = 0; i < points_.size(); i++)
    {
        points_[i].Dump(ofs);
        ofs  << std::endl;
    }

    for (auto iter = order_.begin(); iter != order_.end(); iter++)
    {
        ofs << *iter << std::endl;
    }

    return ofs;
}


void Boundary::GetArrayData(double** x, double **y, double **dx, double **dy, double **ddx, double **ddy) const
{
    *x = (double*) malloc(size()*sizeof(double));
    *y = (double*) malloc(size()*sizeof(double));
    *dx = (double*) malloc(size()*sizeof(double));
    *dy = (double*) malloc(size()*sizeof(double));
    *ddx = (double*) malloc(size()*sizeof(double));
    *ddy = (double*) malloc(size()*sizeof(double));

    #pragma omp parallel for
    for (int i = 0; i < size(); i++)
    {
        (*x)[i] = GetOrderedPoint(i).x;
        (*y)[i] = GetOrderedPoint(i).y;
        (*dx)[i] = GetOrderedPoint(i).D.x();
        (*dy)[i] = GetOrderedPoint(i).D.y();
        (*ddx)[i] = GetOrderedPoint(i).DD.x();
        (*ddy)[i] = GetOrderedPoint(i).DD.y();
    }
}


void Boundary::GeneratePoints(const std::vector<std::vector<bool> >& points, double h)
{
    for (int i = 0; i < points.size(); i++)
    {
        for (int j = 0; j < points[i].size(); j++)
        {
            if (IsBoundary(i, j, points))
            {
                Point b;
                b.i = i;
                b.j = j;
                b.x = h*i;
                b.y = h*j;
                
                points_.push_back(b);
            }         
                
        }
    }

    std::sort(points_.begin(), points_.end());
}


void Boundary::GenerateGeometricalData(const std::vector<std::vector<bool> >& points, double h)
{
    double dt = 2.0*3.14159 / order_.size();

    //compute normal by averaging
    for (int k = 0; k < points_.size(); k++)
    {       
        points_[k].normal = ComputeNormal(points_[k].i, points_[k].j, points);
    }

    //compute the first derivative
    for (auto iter = order_.begin(); iter != order_.end(); iter++)
    {
        auto iterl = iter != order_.begin() ? std::prev(iter) : --order_.end();
        auto iterr = iter != --order_.end() ? std::next(iter) : order_.begin();
        double dx = (points_[*iterr].x - points_[*iterl].x)/(2.0*dt);
        double dy = (points_[*iterr].y - points_[*iterl].y)/(2.0*dt);

        points_[*iter].D = Vec2D(dx, dy);
    }

    //compute the second derivative
    for (auto iter = order_.begin(); iter != order_.end(); iter++)
    {
        auto iterl = iter != order_.begin() ? std::prev(iter) : --order_.end();
        auto iterr = iter != --order_.end() ? std::next(iter) : order_.begin();
        double ddx = (points_[*iterr].D.x() - points_[*iterl].D.x())/(2.0*dt);
        double ddy = (points_[*iterr].D.y() - points_[*iterl].D.y())/(2.0*dt);

        points_[*iter].DD = Vec2D(ddx, ddy);
    }
}


Vec2D Boundary::ComputeNormal(int i, int j, const std::vector<std::vector<bool> >& points) const
{
    double normalx = 0.0;
    double normaly = 0.0;
    int counter = 0;
    for (int ii = i-1; ii <= i+1; ii++)
    {
        for (int jj = j-1; jj <= j+1; jj++)
        {
            if (ii >= 0 && ii < nx_ && jj >= 0 && jj < ny_ &&
                !points[ii][jj])
            {
                double di = ii-i;
                double dj = jj-j;

                normalx += di;
                normaly += dj;
                counter++;
            }
        }
    }

    double norm = std::sqrt(normalx*normalx + normaly*normaly);
    if (norm != 0)
    {
        normalx /= norm;
        normaly /= norm;
    }
    return Vec2D(normalx, normaly);
}


bool Boundary::IsBoundary(int i, int j, const std::vector<std::vector<bool> >& points) const
{
    bool is_boundary = false;
    if (points[i][j]) //interior
    {
        if (i != 0 && !points[i-1][j])
            is_boundary = true;
        
        if (i < points.size()-1 && !points[i+1][j])
            is_boundary = true;

        if (j != 0 && !points[i][j-1])
            is_boundary = true;
        
        if (j < points[i].size()-1 && !points[i][j+1])
            is_boundary = true;
    }

    return is_boundary;
}


std::vector<const Boundary::Point*> Boundary::GetNeighbors(const Point& p, const std::vector<std::vector<bool> >& points) const
{
    std::vector<const Point*> neighbors;

    if (p.j > 0 && IsBoundary(p.i, p.j-1, points))
    {
        Point tmp;
        tmp.i = p.i;
        tmp.j = p.j-1;
        neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
    }
    if (p.j < ny_-1 && IsBoundary(p.i, p.j+1, points))
    {
        Point tmp;
        tmp.i = p.i;
        tmp.j = p.j+1;
        neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
    }
    if (p.i > 0)
    {
        if (IsBoundary(p.i-1, p.j, points))
        {
            Point tmp;
            tmp.i = p.i-1;
            tmp.j = p.j;
            neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
        }
        if (p.j > 0 && IsBoundary(p.i-1, p.j-1, points))
        {
            Point tmp;
            tmp.i = p.i-1;
            tmp.j = p.j-1;
            neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
        }
        if (p.j < ny_-1 && IsBoundary(p.i-1, p.j+1, points))
        {
            Point tmp;
            tmp.i = p.i-1;
            tmp.j = p.j+1;
            neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
        }
    }
    if (p.i < nx_-1)
    {
        if (IsBoundary(p.i+1, p.j, points))
        {
            Point tmp;
            tmp.i = p.i+1;
            tmp.j = p.j;
            neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
        }
        if (p.j > 0 && IsBoundary(p.i+1, p.j-1, points))
        {
            Point tmp;
            tmp.i = p.i+1;
            tmp.j = p.j-1;
            neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
        }
        if (p.j < ny_-1 && IsBoundary(p.i+1, p.j+1, points))
        {
            Point tmp;
            tmp.i = p.i+1;
            tmp.j = p.j+1;
            neighbors.push_back(&*(std::lower_bound(points_.begin(), points_.end(), tmp)));
        }
    }

    return neighbors;
}


void Boundary::GenerateOrdering(const std::vector<std::vector<bool> >& points)
{
    order_.push_back(0); //we start at the first point_.
    std::vector<Point> boundary_copy = points_;
    boundary_copy.erase(boundary_copy.begin());
    
    while (!boundary_copy.empty())
    {
        std::vector<int> possibilities; //if there are multiple potential neighbours.
        possibilities.clear();
        for (int ii = points_[order_.back()].i-1; ii <= points_[order_.back()].i+1; ii++)
        {
            for (int jj = points_[order_.back()].j-1; jj <= points_[order_.back()].j+1; jj++)
            {
                if ((ii == points_[order_.back()].i && jj == points_[order_.back()].j) || //skip the own point
                    !(ii >= 0 && ii < points.size() && jj >= 0 && jj < points[ii].size()))
                    continue;

                Point tmp;
                tmp.i = ii;
                tmp.j = jj;

                if (std::binary_search(boundary_copy.begin(), boundary_copy.end(), tmp))
                {
                    auto it_pts = std::lower_bound(points_.begin(), points_.end(), tmp);
                    int pos = std::distance(points_.begin(), it_pts);
                    possibilities.push_back(pos);
                }
            }
        }

        if (possibilities.empty())
            throw std::runtime_error("Ran into an infinite loop when constructing the boundary");

        //i choose the element that has a normal closest to the normal of the current point
        Vec2D n = ComputeNormal(points_[order_.back()].i, points_[order_.back()].j, points);
        int iter_choice = 0;
        double maximum = -1000000.0;
        for (int k = 0; k < possibilities.size(); k++)
        {
            Vec2D new_n = ComputeNormal(points_[possibilities[k]].i, points_[possibilities[k]].j, points);
            if (maximum < Vec2D::Dot(n, new_n))
            {
                iter_choice = k;
                maximum = Vec2D::Dot(n, new_n);
            }
        }

        auto iter = std::lower_bound(boundary_copy.begin(), boundary_copy.end(), points_[possibilities[iter_choice]]);
        order_.push_back(possibilities[iter_choice]);
        //std::cout << "Node = " << possibilities[iter_choice] << " appended" << std::endl;
        boundary_copy.erase(iter);
    }

    if (order_.size() != points_.size())
        throw std::runtime_error("Not all elements occur in the boundary ordering!");
}