#ifndef VEC2D_H
#define VEC2D_H

class Vec2D
{
public:
    Vec2D(double x, double y);

    double x() const {return x_;}
    double y() const {return y_;}

private:
    double x_;
    double y_;
};


#endif /* VEC2D_H */