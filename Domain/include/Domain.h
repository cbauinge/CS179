#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <tuple>
#include <fstream>

#include "Vec2D.h"

class Domain
{
    friend class DomainTest;
    
public:
////TODO: implement defult constructor
    Domain(std::vector<std::vector<bool> > points);

    const std::vector<std::vector<bool> >& GetPoints() const {return points_;}
    const std::vector<std::pair<int, int> >& GetInterior() const {return interior_;}
    const std::vector<int>& GetBoundary() const {return boundary_;}
    const std::vector<Vec2D>& GetNormals() const {return normals_;}
    double GetH() const {return h_;}

    std::ostream& Dump(std::ostream& ofs) const;


private:
    ///Takes the points stored in points_ and extracts the interior points.
    ///I.e. the points marked 'true' and returns the result.
    std::vector<std::pair<int, int> > GenerateInterior() const;

    ///Generates the boundary using the information stored in points_ AND 
    ///interior_.
    std::vector<int> GenerateBoundary() const;

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

    ///Avector of indices pointing to the interior points which belong to the 
    ///boundary.
    std::vector<int> boundary_;

    ///A vector where the position in the vector is the boundary pixel and the 
    ///value is the normal
    std::vector<Vec2D> normals_;

    ///Discretization size assuming that the whole bmp in x direction covers the interval [0, 1]
    double h_; 
};



#endif /* DOMAIN_H */