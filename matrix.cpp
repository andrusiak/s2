#include <iostream>
#include <iomanip>
#include <cmath>
#include "matrix.h"

using std::ostream;
using std::setprecision;
using std::cout;
using std::endl;

ostream& operator<<(ostream& os, const Matrix& A)
{
    bool useLabel = A.label != "";
    for (unsigned i = 0; i < A.n; i++)
    {
        if (useLabel && i == A.n / 2)
            cout << A.label << " = ";
        os << (useLabel ? "\t" : "") << "[";
        for (unsigned j = 0; j < A.m; j++)
            os << setprecision(3) << A.Mat[i][j] << (j < A.m-1 ? "\t" : "");
        os << "]" << endl;
    }
    return os;
}

Matrix Matrix::operator[](unsigned ind) const
{
    if (this->n <= ind)
        throw Matrix::OperationUndefined("Row index larger that the number of rows");

    Matrix res(1, this->m, "r"+ind);
    for (unsigned i = 0; i < this->m; i++)
        res.Mat[0][i] = this->Mat[ind][i];
    return res;
}

Matrix Matrix::operator+(const Matrix& B) const
{
    if (B.n != n || B.m != m)
        throw Matrix::OperationUndefined("A and B have different sizes");

    Matrix res(n, m, label != "" ? label + "+" + B.label : "");
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < m; j++)
            res.Mat[i][j] = this->Mat[i][j] + B.Mat[i][j];
    return res;
}

Matrix Matrix::operator-(const Matrix& B) const
{
    if (B.n != n || B.m != m)
        throw Matrix::OperationUndefined("A and B have different sizes");

    Matrix res(n, m, label != "" ? label + "-" + B.label : "");
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < m; j++)
            res.Mat[i][j] = this->Mat[i][j] - B.Mat[i][j];
    return res;
}

Matrix& Matrix::operator*=(const double d)
{
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < m; j++)
            this->Mat[i][j] *= d;
    return *this;
}

Matrix Matrix::operator*(const Matrix& B) const
{
    if (m != B.n)
        throw OperationUndefined("A's cols must be same as B's rows");

    Matrix res(n, B.m, label != "" ? label + "*" + B.label : "");
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < B.m; j++)
            for (unsigned k = 0; k < m; k++)
                res(i, j) += this->Mat[i][k] * B.Mat[k][j];
    return res;
}

Matrix Matrix::transpose() const
{
    Matrix res(m, n, label + "^T");
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < m; j++)
            res.Mat[j][i] = this->Mat[i][j];
    return res;
}

Matrix Matrix::inverse() const
{
    if (n != m)
        throw Matrix::OperationUndefined("A not square");

    const double eps = 1.0e-7;
    unsigned n2 = 2*n;
    Matrix temp(n, n2,"T");
    for (unsigned i = 0; i < n; i++)
    {
        for (unsigned j = 0; j < n; j++)
            temp.Mat[i][j] = this->Mat[i][j];
        temp.Mat[i][n+i] = 1;
    }

    for (unsigned i = 0; i < n; i++)
    {
        unsigned maxrow = i;
        for (unsigned j = i+1; j < n; j++)
            if (fabs(temp.Mat[j][i]) > fabs(temp.Mat[maxrow][i]))
                maxrow = j;

        if (fabs(temp.Mat[maxrow][i]) <= eps)
            throw OperationUndefined("A is singular");
        else if (i != maxrow)
            swapRows(temp, i, maxrow);

        double normal = temp.Mat[i][i];
        for (unsigned j = i; j < n2; j++)
            temp.Mat[i][j] /= normal;

        for (unsigned j = i+1; j < n; j++)
        {
            double first = temp.Mat[j][i];
            for (unsigned k = i; k < n2; k++)
                temp.Mat[j][k] -= first * temp.Mat[i][k];
        }

        //cout << temp;
    }

    for (unsigned i = 1; i < n; i++)
        for (unsigned j = 0; j < i; j++)
        {
            double up = temp.Mat[j][i];
            for (unsigned k = i; k < n2; k++)
                temp.Mat[j][k] -= up * temp.Mat[i][k];
        }

    //cout << temp;

    Matrix res(n, n, label != "" ? label + "-1" : "");
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = n; j < n2; j++)
            res.Mat[i][j-n] = fabs(temp.Mat[i][j]) <= eps ? 0.0 : temp.Mat[i][j];
    return res;
}
