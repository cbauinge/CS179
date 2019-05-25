#ifndef VEC2D_H
#define VEC2D_H

#include <fstream>

class Vec2D
{
    friend class Vec2DTest;

public:
    Vec2D(double x = 0, double y = 0);

    double x() const {return x_;}
    double y() const {return y_;}

    Vec2D operator+(const Vec2D& rhs) const;
    Vec2D& operator+=(const Vec2D& rhs);
    Vec2D& operator=(const Vec2D& rhs);

    std::ostream& Dump(std::ostream& ofs) const;


    static double Norm(const Vec2D& v);
    static Vec2D Normalize(const Vec2D& v);

    

private:
    double x_;
    double y_;
};


#endif /* VEC2D_H */