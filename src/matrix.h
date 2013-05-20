#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <string>

using std::ostream;
using std::string;

class Matrix
{
    public:
        struct BadIndex
        {
            string msg;
            BadIndex(string message) : msg(message) {} 
        };

        struct OperationUndefined
        {
            string msg;
            OperationUndefined(string message) : msg(message) {} 
        };

    static const unsigned NMAX = 128;

    private:
        string label;
        double Mat[NMAX][NMAX];
        unsigned n, m;

        static void swapRows(Matrix& mat, unsigned i, unsigned j)
        {
            for (unsigned k = 0; k < mat.m; k++)
            {
                double t = mat.Mat[i][k];
                mat.Mat[i][k] = mat.Mat[j][k];
                mat.Mat[j][k] = t;
            }
        }

    public:
        Matrix(unsigned nn, unsigned mm, const string lbl = "") 
        {
            if (nn > NMAX)
                nn = NMAX;
            if (mm > NMAX)
                mm = NMAX;
            n = nn;
            m = mm;
            label = lbl.substr(0, 4);
            for (unsigned i = 0; i < n; i++)
                for (unsigned j = 0; j < m; j++)
                    this->Mat[i][j] = 0.0;
        }

        inline string getLabel() { return label; }
        inline unsigned getRows() { return n; }
        inline unsigned getCols() { return m; }
        inline void setLabel(string l) { label = l; }

        inline double& operator() (unsigned row, unsigned col)
        {
            if (row < n && col < m)
                return Mat[row][col];
            throw Matrix::BadIndex("index out of bounds");
        }

        inline double operator() (unsigned row, unsigned col) const
        {
            if (row < n && col < m)
                return Mat[row][col];
            throw Matrix::BadIndex("out of bounds");
        }

        Matrix& operator*=(const double d);
        Matrix operator+(const Matrix& B) const;
        Matrix operator-(const Matrix& B) const;
        Matrix operator*(const Matrix& B) const;
        Matrix operator[](unsigned ind) const;
        Matrix transpose() const;
        Matrix inverse() const;

        friend ostream& operator<<(ostream& os, const Matrix& mat);
};


#endif

