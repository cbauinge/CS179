#include "SolverEigen.h"

#include <cmath>
#include <iostream>


using namespace std::complex_literals;

std::vector<double> SolverEigen::Solve(const Domain& dom) const
{
    Matrix A = SetupA(dom);
    Matrix I = Matrix::Identity(A.rows(), A.cols());
    //std::cout << "A = " << A << std::endl;
    Vector b = SetupRhs(dom);

    //solve for density
    Eigen::ColPivHouseholderQR<Matrix> dec(I - A);
    Vector x = dec.solve(b);
    //std::cout << std::endl << "x = " << x << std::endl;

    //now, evaluate the result on all interior points. For that, setup the evaluation matrix
    Matrix Eval = SetupEval(dom);
    //std::cout << "Eval = " << Eval << std::endl;

    //Evaluate the result in the inteiror (wo bounadary)
    Vector evaluation = Eval*x;
    //std::cout << "FInal inner result = " << evaluation << std::endl;
    
    std::vector<double> result;
    result.resize(evaluation.size());
    for (int i = 0; i < result.size(); i++)
    {
        result[i] = evaluation(i).real();
    }

    return result;
}


Matrix SolverEigen::SetupA(const Domain& dom) const
{
    int n = dom.GetBoundary().size();
    Matrix A(n, n);
    double dt = 2.0*M_PI/n;
    int fn = floor(n/2);


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A(i, j) = ComputeR(dt*j, dt*i, fn)*K1(i, j, dom) + dt*K2(i, j, dom);
        }
    }


    return A;
}


Vector SolverEigen::SetupRhs(const Domain& dom) const
{
    int n = dom.GetBoundary().size();
    Vector b(n);

    for (int i = 0; i < n; i++)
    {
        b(i) = -2.0*bc_[i];
    }

    return b;
}


Matrix SolverEigen::SetupEval(const Domain& dom) const
{
    Matrix Eval(dom.GetInteriorWOBoundary().size(), dom.GetBoundary().size());
    
    for (int i = 0; i < Eval.rows(); i++)
    {
        for (int j = 0; j < Eval.cols(); j++)
        {   
            std::pair<int, int> coordinates = dom.GetInteriorWOBoundary()[i];
            Eval(i, j) = Lxy(coordinates.first*dom.GetH(), coordinates.second*dom.GetH(), j, dom);
        }
    }

    return Eval;
}
    

Matrix SolverEigen::SetupR(const Domain& dom) const
{
    int n = dom.GetBoundary().size();
    Matrix R(n, n);

    int fn = std::floor(n/2);

    //double dt = 2*PI/n;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            R(i, j) = 0.0;
            for (int k = 1; k < fn; k++)
            {
                R(i, j) += 1.0/k*cos(k*i*M_PI/double(fn)) + std::pow(-1.0, i)/(2.0*fn);
            }
            R(i, j) /= -fn;
        }
    }

    return R;    
}

double SolverEigen::ComputeR(double tj, double t, int n) const
{
    double sum = 0.0;

    for (int m = 1; m < n; m++)
    {
        sum += 1.0/m * cos(m*(t - tj)) + 1/(2.0*n)*cos(n*(t - tj)); 
    }
    sum /= -n;
}


complex SolverEigen::K(int i, int j, const Domain& dom) const
{
    //double dt = 2.0*M_PI/dom.GetBoundary().size();

    complex val = 0.0;

    if (abs(i-j) < 0.5)
    {
        Vec2D D = dom.GetBoundary().GetOrderedPoint(i).D;
        double dx = D.x();
        double dy = D.y();
        double ddx = dom.GetBoundary().GetOrderedPoint(i).DD.x();
        double ddy = dom.GetBoundary().GetOrderedPoint(i).DD.y();
        val = -1.0/(2.0*M_PI*Vec2D::Dot(D, D))*(dx*ddy - dy*ddx);
    }
    else
    {
        Vec2D xji(dom.GetBoundary().GetOrderedPoint(j).x - dom.GetBoundary().GetOrderedPoint(i).x,
            dom.GetBoundary().GetOrderedPoint(j).y - dom.GetBoundary().GetOrderedPoint(i).y);
        double normxji = Vec2D::Norm(xji);
        Vec2D D = dom.GetBoundary().GetOrderedPoint(j).D;
        val = -complex(1i*k_/(2.0*normxji)*(D.y()*xji.x() - D.x()*xji.y()))*Hankel(1, k_*normxji);
    }

    return val;
}


complex SolverEigen::K1(int i, int j, const Domain& dom) const
{
    complex val = 0.0;

    if (abs(i-j) > 0.5)
    {
        Vec2D xij(dom.GetBoundary().GetOrderedPoint(i).x - dom.GetBoundary().GetOrderedPoint(j).x,
                dom.GetBoundary().GetOrderedPoint(i).y - dom.GetBoundary().GetOrderedPoint(j).y);
        double normxij = Vec2D::Norm(xij);
        Vec2D D = dom.GetBoundary().GetOrderedPoint(j).D;
        val =  k_/(2.0*M_PI*normxij)*(D.y()*xij.x() - D.x()*xij.y())*std::cyl_bessel_j(1, k_*normxij);
    }
    return val;
}


complex SolverEigen::K2(int i, int j, const Domain& dom) const
{
    complex val = 0.0;
    if (abs(i-j) > 0.5)
    {
        double dt = 2.0*M_PI/dom.GetBoundary().size();
        double tmp = sin((i*dt - j*dt)/2.0);
        val = K(i, j, dom) - K1(i, j, dom)*log(4*tmp*tmp);
    }
    else
    {
        val = K(i, j, dom);
    }

    return val;
}


complex SolverEigen::Lxy(double x, double y, int i, const Domain& dom) const
{
    Vec2D dist(x - dom.GetBoundary().GetOrderedPoint(i).x, y - dom.GetBoundary().GetOrderedPoint(i).y);
    double normdist = Vec2D::Norm(dist);
    Vec2D D = dom.GetBoundary().GetOrderedPoint(i).D;
    double dt = 2.0*M_PI/dom.GetBoundary().size();
    complex val = complex(1i * k_ * dt /(4.0*normdist)*(D.y()*dist.x() - D.x()*dist.y()))*Hankel(1, k_*normdist);
    return val;
}