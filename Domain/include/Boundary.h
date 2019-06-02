#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include <list>
#include <fstream>
#include "Vec2D.h"


/// @brief Class that represents the boundary of a domain.
/// 
/// Consists of boundary points which are ordered according to their position in the
/// domain and an indexing set which gives the position in the boundary.
class Boundary
{
    friend class BoundaryTest;

public:
    /// @brief Point struct representing one point on the boundary.
    struct Point
    {
        int i; ///< index of x position in the grid
        int j; ///< index of y position in the grid
        double x; ///< x position
        double y; ///< y position
        Vec2D normal; ///< normal vector pointing outwards
        Vec2D D; ///< First derivative at this point (in the parametrization)
        Vec2D DD; ///< Second derivativ at this point.

        /// @brief Comparison operator for sorting.
        bool operator<(const Point& rhs) const
        {
            return (i < rhs.i ) || (i == rhs.i && j < rhs.j);
        }

        /// @brief Function that writes the struct to the given ofs. 
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
    /// @brief Constructor which generates the boundary from a givne bitmap.
    /// @param [in] points. Matrix of the bitmap.
    /// @param [in] h. size between cells.
    Boundary(const std::vector<std::vector<bool> >& points, double h);

    /// @brief Getter for the boundayr points.
    const std::vector<Point>& GetPoints() const {return points_;}

    /// @brief Getter for the number of boundary points.
    int size() const {return points_.size();}

    /// @brief Function that gets point i in the boundary.
    /// Getting points i-1 and i+1 too ensures that thei are nighbours in the boundary.
    const Point& GetOrderedPoint(int i) const {return points_[order_[i]];}

    /// @brief Write the boundary to the given ostream.
    /// Used for visuaization.
    std::ostream& Dump(std::ostream& ofs) const;

    /// @brief Funciton that check if the point at (i, j) is a boundary point.
    /// NOTE this is bad design
    bool IsBoundary(int i, int j, const std::vector<std::vector<bool> >& points) const;

    ///@brief function that writes the points in given arrays.
    void GetArrayData(double** x, double **y, double **dx, double **dy, double **ddx, double **ddy) const;
private:
    /// @brief Generates the boundary points (i, j, x, y) from the given bitmap and distance h.
    void GeneratePoints(const std::vector<std::vector<bool> >& points, double h);

    /// @brief Function that computes the geometrical data, normal, D, and DD.
    ///
    /// Function that assumes that boundayr points are setup with information i, j, x, y, 
    /// stored in points_ and sorted by operator<. Then computes D and DD (first and second derivatives)
    /// for each point 
    void GenerateGeometricalData(const std::vector<std::vector<bool> >& points, double h);

    /// @brief FuncitonThat computes the normal of a given point i, j according to the bitmap points. 
    Vec2D ComputeNormal(int i, int j, const std::vector<std::vector<bool> >& points) const;

    /// @brief Get all the neighbouring boundary points.
    std::vector<const Boundary::Point*> GetNeighbors(const Point& p, const std::vector<std::vector<bool> >& points) const;

    /// @brief Generate the order of the particles in th eboundary.
    void GenerateOrdering(const std::vector<std::vector<bool> >& points);

private:
    /// all the points in the boundary. They should always be sorted by their position.
    std::vector<Point> points_;

    ///the ordering. neighbors in this vector ar eneighbors on the boundayr. the value is the position in points_;
    std::vector<int> order_; 
    int nx_; ///<number of elements in x direction in the bitmap.
    int ny_; ///<number of elements in y direction in the bitmap.
};

#endif /* BOUNDARY_H */