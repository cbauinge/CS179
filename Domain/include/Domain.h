#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <tuple>

#include "Vec2D.h"

class Domain
{
    friend class DomainTest;
    
public:
    Domain(std::vector<std::vector<bool> > points);

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
};


#endif /* DOMAIN_H */