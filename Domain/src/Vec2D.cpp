#include "Vec2D.h"

#include<cmath>
#include<iostream>


Vec2D::Vec2D(double x, double y) :
    x_(x), y_(y)
{
    
}


Vec2D Vec2D::operator+(const Vec2D& rhs) const
{
    Vec2D val(rhs.x() + x_, rhs.y() + y_);
    return val;
}


Vec2D& Vec2D::operator+=(const Vec2D& rhs)
{
    *this = *this + rhs;
    return *this;
}


Vec2D& Vec2D::operator=(const Vec2D& rhs)
{
    if (this == &rhs)
        return *this;

    x_ = rhs.x();
    y_ = rhs.y();

    return *this;
}


std::ostream& Vec2D::Dump(std::ostream& ofs) const
{
    ofs << x() << ", " << y() << std::endl;

    return ofs;
}


double Vec2D::Norm(const Vec2D& v)
{
    return std::sqrt(v.x()*v.x() + v.y()*v.y());
}


Vec2D Vec2D::Normalize(const Vec2D& v)
{
    double norm = Norm(v);
    Vec2D val(0, 0);
    if (norm != 0.0)
        val = Vec2D(v.x()/norm, v.y()/norm);
    else
        std::cout << "Warning: vector to normalize is (0,0)" << std::endl;
    
    return val;
}