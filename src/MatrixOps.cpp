#include <stdexcept>
#include <algorithm> //for std::find
#include <cmath>
#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/MatrixOps.h"

namespace la {

void PrintMatrix(const Matrix& m) {
    for (const Vector& v : m) {
        PrintVector(v);
        std::cout << "\n";
    }
}

Matrix MatrixAddition(const Matrix& m1, const Matrix& m2) {
    std::vector<Matrix::size_type> m1Dim = m1.dimension();
    std::vector<Matrix::size_type> m2Dim = m2.dimension();

    if (m1Dim != m2Dim) {
        throw std::invalid_argument("Matrices must have same dimension.");
    } 

    Matrix ma = Matrix(m1Dim[0],m1Dim[1]); // Initialize zero matrix we will modify for final output

    for (int i = 0; i < m1Dim[0]; ++i) {
        for (int j = 0; j < m1Dim[1]; ++j) {
            ma.row(i)[j] = m1.row(i)[j] + m2.row(i)[j];
        }
    }

    return ma;
}

Matrix MatrixSubtraction(const Matrix& m1, const Matrix& m2) {
    Matrix m2Neg = ScalarMultiplication(m2,-1);
    return MatrixAddition(m1,m2Neg);
}

Matrix ScalarMultiplication(const Matrix& m, double k) {
    Matrix sm = Matrix(m);

    for (Vector& r : sm) {
        for (double& c : r) {
            c *= k;
        }
    }

    return sm;
}

// Matrix MatrixMultiplication(const Matrix& m1, const Matrix& m2) {
//     std::vector<Matrix::size_type> m1Dim = m1.dimension();
//     std::vector<Matrix::size_type> m2Dim = m2.dimension();

//     if (m1Dim[1] != m2Dim[0]) {
//         throw std::invalid_argument("m1 column dimension must equal m2 row dimension.");
//     }

//     Matrix mm = Matrix(m1Dim[0],m2Dim[1]);

//     for (Matrix::size_type i = 0; i < m1Dim[0]; ++i) {
//         Vector r = Vector(m2Dim[1],0.0);
//         for (Matrix::size_type j = 0; j < m2Dim[1]; ++j) {
//             double entry = 0.0;
//             for (Matrix::size_type k = 0; k < m1Dim[1]; ++k) {
//                 entry += m1.row(i)[k] * m2.row(k)[j];
//             }
//             r[j] = entry;
//         }
//         mm.row(i) = r;
//     }

//     return mm;
    
// }

// double Determinant2x2(const Matrix& m) {
//     std::vector<Matrix::size_type> mDim = m.dimension();
//     if (mDim[0] != 2 || mDim[1] != 0) {
//         throw new std::invalid_argument("Must enter a 2x2 matrix.");
//     }

//     return (m.row(0)[0] * m.row(1)[1]) - (m.row(0)[1] * m.row(1)[0]);
// }

// Matrix MinorMatrix(const Matrix& m, const std::set<Matrix::size_type>& r, const std::set<Matrix::size_type>& c) {
//     /*
//     Parameters:
//     -----------
//     -> m: -> The matrix we wish to find Minor Matrix for
//     -> r -> index(s) for the delete row
//     -> c -> index(s) for the delete column

//     Note: A MinorMatrix is a smaller matrix obtained from a larget one by deleting one or more rows and/or one columns
//     */

//     std::vector<Matrix::size_type> mDim = m.dimension();
//     std::size_t rSize = r.size();
//     std::size_t cSize = c.size();
//     Matrix::size_type mmRow = mDim[0] - rSize;
//     Matrix::size_type mmCol = mDim[1] - cSize;

//     Matrix mm = Matrix(mmRow,mmCol); // Initialize zero matrix we will fill in

//     Matrix::size_type curr_mmRow = 0;
//     for (Matrix::size_type i = 0; i < mDim[0]; ++i) {
//         if (r.find(i) != r.end()) {
//             continue;
//         }
//         Matrix::size_type curr_mmCol = 0;
//         for (Matrix::size_type j = 0; j < mDim[1]; ++j) {
//             if (c.find(j) != c.end()) {
//                 continue;
//             }
//             mm.row(curr_mmRow)[curr_mmCol] = m.row(i)[j];
//             ++curr_mmCol;
//         }
//         ++curr_mmRow;
//     }

//     return mm;

// }

// double Determinant(const Matrix& m) {
//     std::vector<Matrix::size_type> mDim = m.dimension();
//     if (mDim[0] != mDim[1]) {
//         throw std::invalid_argument("Matrix must be a square matrix to compute the determinant.");
//     }

//     if (mDim[0] == 1) {
//         return m.row(0)[0]; // 1x1 determinant has matrix equal to its only entry
//     }
//     else if (mDim[0] == 2) {
//         return Determinant2x2(m);
//     }
//     else { // Use Cofactor Expansion
//         Matrix::size_type r_delete = 0;             // row to delete for MinorMatrix
//         std::vector<Matrix::size_type> z_indices;   // indicates col index in cofactor row where entry = 0
//         Matrix::size_type z_count = 0;                            // indicates # of zeros in cofactor row

//         // Below for loop searches for row with the most zeros to use as delete row in MinorMatrix
//         for (Matrix::size_type i = 0; i < mDim[0]; ++i) {
//             if (z_count == mDim[1]) {
//                 return 0.0; // If zero matrix, determinant = 0
//             }

//             Matrix::size_type r_num_zero = 0; // Current # of 0 entries in row
//             std::vector<Matrix::size_type> r_z_indices; // col index for current row check 0 entries

//             for (Matrix::size_type j = 0; j < mDim[1]; ++j) {
//                 if (m.row(i)[j] == 0) {
//                     ++r_num_zero;
//                     r_z_indices.push_back(j);
//                 }
//             }

//             if (r_num_zero > z_count) {
//                 r_delete = i;
//                 z_count = r_num_zero;
//                 z_indices = std::move(r_z_indices);
//             }
//         }

//         // Now we search for the first nonzero entry in our row to serve as delete col in MinorMatrix
//         double det = 0.0;
//         for (Matrix::size_type k = 0; k < mDim[1]; ++k) {
//             if (std::find(z_indices.begin(), z_indices.end(), k) != z_indices.end()) {
//                 continue;
//             }
//             else {
//                 Matrix mm = MinorMatrix(m, std::set<Matrix::size_type>{r_delete}, std::set<Matrix::size_type>{k});

//                 double sign = ((r_delete + k) & 1) ? -1.0 : 1.0;    // bitwise parity
//                 det += std::pow(-1,(r_delete + k)) * m.row(r_delete)[k] * Determinant(mm);
//             }
//         }

//         return det;

//     }
    
// }

// Matrix RowScale(const Matrix& m, Matrix::size_type& r) {
//     double f = 0.0;
//     for (double x : m.row(r)) {
//         if (x != 0) {
//             f = x;
//             break;
//         }
//     }

//     if (f == 0.0) {
//         return m;
//     }
//     else {
//         Matrix m2(m);
    
//         for (double& c : m2.row(r)) {
//             c /= f;
//         }
//     }
// }

// Matrix RowSwap(const Matrix& m, Matrix::size_type& r1, Matrix::size_type& r2) {
//     Matrix rs(m);
//     return std::swap(rs.row(r1),rs.row(r2));
// }


// Matrix RowReplacement(const Matrix&m, Matrix::size_type& r1, Matrix::size_type& r2);
// Matrix RowEchelonForm(const Matrix& m);
// Matrix ReducedRowEchelonForm(const Matrix& m);

}