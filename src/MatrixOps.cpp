#include <stdexcept>
#include <iostream>
#include "claire/la/MatrixOps.h"

namespace la {

Matrix MatrixMultiplication(const Matrix& m1, const Matrix& m2) {
    std::vector<Matrix::size_type> m1Dim = m1.dimension();
    std::vector<Matrix::size_type> m2Dim = m2.dimension();

    if (m1Dim[1] != m2Dim[0]) {
        throw std::invalid_argument("m1 column dimension must equal m2 row dimension.");
    }

    Matrix mm = Matrix(m1Dim[0],m2Dim[1]);

    for (Matrix::size_type i = 0; i < m1Dim[0]; ++i) {
        Vector r = Vector(m2Dim[1],0.0);
        for (Matrix::size_type j = 0; j < m2Dim[1]; ++j) {
            double entry = 0.0;
            for (Matrix::size_type k = 0; k < m1Dim[1]; ++k) {
                entry += m1.row(i)[k] * m2.row(k)[j];
            }
            r[j] = entry;
        }
        mm.row(i) = r;
    }

    return mm;
    
}

}