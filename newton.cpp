#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "matrix.h"

using namespace std;

unsigned M = 99, N = 2;

double H(const Matrix& theta, const Matrix& xs)
{
    Matrix t = theta * xs;
    return t(0,0);
}

double logistic(double h)
{
    return 1.0 / (1.0 + exp(-h));
}

void newton(Matrix& theta, const Matrix& X, const Matrix& y)
{
    for (int kk = 0; kk < 10; kk++)
    {
        Matrix deriv(N+1, 1, "der");
        for (unsigned j = 0; j < N+1; j++)
        {
            double sum = 0.0;
            for (unsigned m = 0; m < M; m++)
                sum += (y(m, 0) - logistic(H(theta.transpose(), X[m].transpose()))) * X[m](0, j);
            deriv(j, 0) = sum;
        }

        double sum = 0.0;
        for (unsigned m = 0; m < M; m++)
        {
            double g = logistic(H(theta.transpose(), X[m].transpose()));
            sum += g*g - g;
        }
        Matrix hessian(N+1, N+1, "H");
        Matrix XtX = X.transpose() * X;
        for (unsigned j = 0; j < N+1; j++)
            for (unsigned k = 0; k < N+1; k++)
                hessian(j, k) = XtX(j, k) * sum;

        theta = theta - (hessian.inverse() * deriv);
        cout << theta << endl;
    }
}

int lllmain()
{
    ifstream q1x("q1x.dat");
    ifstream q1y("q1y.dat");
    if (!q1x.is_open() || !q1y.is_open())
    {
        cerr << "Failed to open file" << endl;
        return 1;
    }


    Matrix X(M, N+1, "X");
    Matrix y(M, 1, "y");
    for (unsigned i = 0; i < M; i++)
    {
        X(i, 0) = 1;
        for (unsigned j = 1; j <= N; j++)
            q1x >> X(i, j);
        q1y >> y(i, 0);
    }
    q1x.close();
    q1y.close();

    Matrix theta(N+1, 1, "Th");
    newton(theta, X, y);

    //Matrix Xt = X.transpose();
    //Matrix theta = (Xt * X).inverse() * Xt * y;
    //theta.setLabel("theta");
    //cout << theta << endl;

    q1x.open("q1x.dat");
    for (unsigned mm = 1; mm <= M; mm++)
    {
        Matrix xs(N+1, 1, "x");
        xs(0,0) = 1;
        //cout << "Input " << N << " features:" << endl;
        for (unsigned i = 1; i <= N; i++)
            q1x >> xs(i, 0);

        double res = logistic(H(theta.transpose(), xs));
        if ((mm < 51 && res > 0.5) || (mm >= 51 && res < 0.5))
            cout << "Logistic(m=" << mm << ") = " << (res > 0.5 ? 1.0 : 0.0) << " -> " << res << endl;
    }

    return 0;
}

