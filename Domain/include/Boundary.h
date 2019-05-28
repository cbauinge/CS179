#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include <list>
#include <fstream>
#include "Vec2D.h"

class Boundary
{
    friend class BoundaryTest;

public:
    struct Point
    {
        int i, j;
        double x, y;
        Vec2D normal;
        Vec2D D; //first derivative
        Vec2D DD; //second derivative

        bool operator<(const Point& rhs) const
        {
            return (i < rhs.i ) || (i == rhs.i && j < rhs.j);
        }

        std::ostream& Dump(std::ostream& ofs) const 
        {
            ofs << "(i, j) = (" << i << ", " << j << ")" << std::endl;
            ofs << "(x, y) = (" << x << ", " << y << ")" << std::endl;
            D.Dump(ofs);
            DD.Dump(ofs);
            return ofs;
        }
    };

public:
    Boundary(const std::vector<std::vector<bool> >& points, double h);

    const std::vector<Point>& GetPoints() const {return points_;}

    int size() const {return points_.size();}
    const Point& operator[](int i) const {return points_[i];}
    const Point& GetOrderedPoint(int i) const {return points_[order_[i]];}
    std::ostream& Dump(std::ostream& ofs) const;


    //funciton that check if the point at (i, j) is a boundary point.
    //this is very bad design
    bool IsBoundary(int i, int j, const std::vector<std::vector<bool> >& points) const;



private:
    void GeneratePoints(const std::vector<std::vector<bool> >& points, double h);

    //function that assumes that boundayr points are setup with information i, j, x, y, 
    //stored in points_ and sorted by operator<. Then computes D and DD (first and second derivatives)
    //for each point 
    void GenerateGeometricalData(const std::vector<std::vector<bool> >& points, double h);

    Vec2D ComputeNormal(int i, int j, const std::vector<std::vector<bool> >& points) const;

    //Get all the neighbouring boundary points
    std::vector<const Boundary::Point*> GetNeighbors(const Point& p, const std::vector<std::vector<bool> >& points) const;

    void GenerateOrdering(const std::vector<std::vector<bool> >& points);

private:
    //all the points in the boundary. They should always be sorted by their position.
    std::vector<Point> points_;
    //the ordering. neighbors in this vector ar eneighbors on the boundayr. the value is the position in points_;
    std::vector<int> order_; 
    int nx_;
    int ny_;
};

#endif /* BOUNDARY_H */