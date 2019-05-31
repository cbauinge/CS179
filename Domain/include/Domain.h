#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <tuple>
#include <fstream>

#include "Vec2D.h"
#include "Boundary.h"


///@brief class representing the domain on which to solve the PDE.
class Domain
{
    friend class DomainTest;
    
public:
////TODO: implement defult constructor

    ///@brief Constructor.
    /// Takes a bitmap and generates the interior points, the boundary points, etc. from it.
    Domain(std::vector<std::vector<bool> > points);

    ///@brief Getter for all the points (inside and outside)
    const std::vector<std::vector<bool> >& GetPoints() const {return points_;}

    ///@brief Getter for the Interior points.
    const std::vector<std::pair<int, int> >& GetInterior() const {return interior_;}

    ///@brief Getter for the interior points excluding boundary points.
    const std::vector<std::pair<int, int> >& GetInteriorWOBoundary() const {return interior_without_boundary_;}

    ///@brief Gettter for the boundary.
    const Boundary& GetBoundary() const {return boundary_;}

    ///@brief Getter for the discretization size.
    double GetH() const {return h_;}

    ///@brief Writes the domain to the given ostream in text format.
    std::ostream& Dump(std::ostream& ofs) const;


private:
    ///Takes the points stored in points_ and extracts the interior points.
    ///I.e. the points marked 'true' and returns the result.
    std::vector<std::pair<int, int> > GenerateInterior() const;

    ///@brief Generate interior points without the boundary.
    std::vector<std::pair<int, int> > GenerateInteriorWOBoundary() const;

    ///Generates the normals using the information stored in points_, interior_
    ///AND boundary_.
    std::vector<Vec2D> GenerateNormals() const;
    
private:
    ///A vector of vectors of points describing the bitmap
    std::vector<std::vector<bool> > points_; 

    ///A vector of pairs describing the pixels belonging to the domain. 
    ///The pair are x and y coordinate in pixels (correspong to points_).
    ///The interior includes the boundary
    std::vector<std::pair<int, int> > interior_;

    ///interior_ without the boundayr elements
    std::vector<std::pair<int, int> > interior_without_boundary_;

    ///Discretization size assuming that the whole bmp in x direction covers the interval [0, 1]
    double h_;
    
    /// boundary
    Boundary boundary_;
};



#endif /* DOMAIN_H */