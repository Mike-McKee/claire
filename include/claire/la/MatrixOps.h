#pragma once
#include <string>
#include "Matrix.h"

namespace la {

    Matrix MatrixAddition(const Matrix& m1, const Matrix& m2);
    Matrix MatrixSubtraction(const Matrix& m1, const Matrix& m2);
    Matrix ScalarMultiplication(const Matrix& m, double k);
    Matrix MatrixMultiplication(const Matrix& m1, const Matrix& m2);
    double Determinant2x2(const Matrix& m);
    Matrix MinorMatrix(const Matrix& m, Matrix::size_type r, Matrix::size_type c);
    double Determinant(const Matrix& m);
    Matrix::size_type Rank(const Matrix& m);
    bool isFullRank(const Matrix& m);
    Matrix InverseMatrix(const Matrix & m);
    Matrix Transpose(const Matrix& m);
    bool isOrthogonal(const Matrix& m1, const Matrix& m2);
    Vector Diagonal(const Matrix& m);
    double Trace(const Matrix& m);
    bool isSymmetric(const Matrix& m);
    bool isDiagonal(const Matrix& m);
    bool isTriangular(const Matrix& m);
    std::string TriangularType(const Matrix& m);

}