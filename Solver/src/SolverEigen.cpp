#include "SolverEigen.h"

#include <cmath>


using namespace std::complex_literals;

std::vector<double> SolverEigen::Solve(const Domain& dom) const
{
    Matrix A = SetupA(dom);
    Vector b = SetupRhs(dom);

    //solve for density
    Eigen::ColPivHouseholderQR<Matrix> dec(Matrix::Identity(A.rows(), A.cols()) - A);
    Vector x = dec.solve(b);

    //now, evaluate the result on all interior points. For that, setup the evaluation matrix
    Matrix Eval(dom.GetInterior().size() - dom.GetBoundary().size(), A.cols());

    //Evaluate the result in the inteiror (wo bounadary)
    Vector evaluation = Eval*x;
    
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
    Matrix R = SetupR(dom);
    double dt = 2.0*M_PI/dom.GetBoundary().size();


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A(i, j) = R(i, j)*K1(i, j, dom) + dt*K2(i, j, dom);
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
        val = 1.0/(2.0*M_PI*Vec2D::Dot(D, D))*(dx*ddy - dy*ddx);
    }
    else
    {
        Vec2D xji(dom.GetBoundary().GetOrderedPoint(j).x - dom.GetBoundary().GetOrderedPoint(i).x,
            dom.GetBoundary().GetOrderedPoint(j).y - dom.GetBoundary().GetOrderedPoint(i).y);
        double normxji = Vec2D::Norm(xji);
        Vec2D D = dom.GetBoundary().GetOrderedPoint(j).D;
        double dx = D.x();
        double dy = D.y();
        val = complex(1i*k_/(2.0*normxji)*(dy*xji.x() - dx*xji.y()))*Hankel(1, k_*normxji);
    }

    return val;
}


complex SolverEigen::K1(int i, int j, const Domain& dom) const
{
    Vec2D xij(dom.GetBoundary().GetOrderedPoint(i).x - dom.GetBoundary().GetOrderedPoint(j).x,
            dom.GetBoundary().GetOrderedPoint(i).y - dom.GetBoundary().GetOrderedPoint(j).y);
    double normxij = Vec2D::Norm(xij);
    Vec2D D = dom.GetBoundary().GetOrderedPoint(j).D;
    double dx = D.x();
    double dy = D.y();
    return k_/(2.0*M_PI*normxij)*(dy*xij.x() - dx*xij.y())*std::cyl_bessel_j(1, k_*normxij);
}


complex SolverEigen::K2(int i, int j, const Domain& dom) const
{
    complex val = 0.0;
    if (abs(i-j) < 0.5)
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
    Vec2D normal = dom.GetBoundary().GetOrderedPoint(i).normal;
    double dt = 2.0*M_PI/dom.GetBoundary().size();
    complex val = complex(1i * k_/4.0*dt*(Vec2D::Dot(dist, normal)/normdist))*Hankel(1, k_*normdist);
    return val;
}