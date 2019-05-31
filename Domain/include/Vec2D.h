#ifndef VEC2D_H
#define VEC2D_H

#include <fstream>

///@brief class representing a 2D vector.
class Vec2D
{
    friend class Vec2DTest;

public:
    ///@brief constructor taking x and y coordinate.
    Vec2D(double x = 0, double y = 0);

    ///@brief Getter for x coordinate.
    double x() const {return x_;}

    ///@brief Getter for y coordinate.
    double y() const {return y_;}

    ///@brief operator +.
    Vec2D operator+(const Vec2D& rhs) const;

    ///@brief add assignment operator.
    Vec2D& operator+=(const Vec2D& rhs);

    ///@brief copy assignment operator.
    Vec2D& operator=(const Vec2D& rhs);

    ///@brief Write vector to ostream.
    std::ostream& Dump(std::ostream& ofs) const;

    ///@brief static function to compute the 2-norm of a vector.
    static double Norm(const Vec2D& v);

    ///@brief static function computing the dot product of two vectors.
    static double Dot(const Vec2D& v1, const Vec2D& v2); 

    ///@brief static function that normalizes a given vector.
    ///@return A vector with norm 1 but the same direction as the given vector.
    static Vec2D Normalize(const Vec2D& v);

    

private:
    double x_; ///<coordinate x.
    double y_; ///<coordinate y.
};


#endif /* VEC2D_H */