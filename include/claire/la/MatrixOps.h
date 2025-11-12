#pragma once
#include <set>
#include <string>
#include "Matrix.h"

/*
FUTURE TODO:
------------
Consider making functions directly modify the input matrix.
Current method of creating multiple copy-vector objects within functions can become expensive.
*/

namespace la {

    void PrintMatrix(const Matrix& m);
    
    bool isSquareMatrix(const Matrix& m);
    bool isZeroMatrix(const Matrix& m);
    bool isIdentityMatrix(const Matrix& m);

    Matrix MatrixAddition(const Matrix& m1, const Matrix& m2);
    Matrix MatrixSubtraction(const Matrix& m1, const Matrix& m2);
    Matrix ScalarMultiplication(const Matrix& m, double k);
    Matrix MatrixMultiplication(const Matrix& m1, const Matrix& m2);
    double Determinant2x2(const Matrix& m);
    Matrix MinorMatrix(const Matrix& m, const std::set<Matrix::size_type>& r, const std::set<Matrix::size_type>& c); // allow multiple delete rows/columns
    double Determinant(const Matrix& m);
    
    Matrix RowScale(const Matrix& m, const Matrix::size_type& r);
    Matrix RowSwap(const Matrix& m, const Matrix::size_type& r1, const Matrix::size_type& r2);
    Matrix RowReplacement(const Matrix&m, const Matrix::size_type& r1, const Matrix::size_type& r2);
    Matrix RowEchelonForm(const Matrix& m);
    Matrix ReducedRowEchelonForm(const Matrix& m);
    
    Matrix::size_type Rank(const Matrix& m);
    bool isFullRank(const Matrix& m);
    Matrix AugmentedMatrix(const Matrix& m1, const Matrix& m2);
    enum class AmSide { L, R };
    Matrix AugmentedMatrixSplit(const Matrix& m, Matrix::size_type r, AmSide side);
    Matrix InverseMatrix(const Matrix & m);
    Matrix Transpose(const Matrix& m);
    bool isOrthogonal(const Matrix& m1, const Matrix& m2);
    bool isOrthonormal(const Matrix&m);
    Vector Diagonal(const Matrix& m);
    double Trace(const Matrix& m);
    bool isSymmetric(const Matrix& m);
    bool isDiagonal(const Matrix& m);
    bool isTriangular(const Matrix& m);
    enum class TriangularMatrixType { Upper, Lower };
    std::vector<std::string> TriangularType(const Matrix& m);

}