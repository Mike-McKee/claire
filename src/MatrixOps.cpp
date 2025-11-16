#include <stdexcept>
#include <algorithm> //for std::find
#include <cmath>
#include <tuple>
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


bool isSquareMatrix(const Matrix& m) {
    if (m.rowCount() == m.colCount()) {
        return true;
    }
    else {
        return false;
    }
}

bool isZeroMatrix(const Matrix& m) {
    for (const Vector& r : m) {
        if (!isZeroVector(r)) {
            return false;
        }
    }

    return true;
}

bool isIdentityMatrix(const Matrix& m) {
    if (!isSquareMatrix(m)) {
        return false;
    }

    std::vector<Matrix::size_type> mDim = m.dimension();

    for (Matrix::size_type i = 0; i < mDim[0]; ++i) {
        for (Matrix::size_type j = 0; j < mDim[1]; ++j) {
            if (i != j && m[i][j] != 0.0) {
                return false;
            }
            else if (i == j && m[i][j] != 1.0) {
                return false;
            }
        }
    }

    return true;
}

double Determinant2x2(const Matrix& m) {
    std::vector<Matrix::size_type> mDim = m.dimension();
    if (mDim[0] != 2 || mDim[1] != 2) {
        throw std::invalid_argument("Must enter a 2x2 matrix.");
    }

    return (m[0][0] * m[1][1]) - (m[0][1] * m[1][0]);
}

Matrix MinorMatrix(const Matrix& m, const std::set<Matrix::size_type>& r, const std::set<Matrix::size_type>& c) {
    /*
    Parameters:
    -----------
    -> m: -> The matrix we wish to find Minor Matrix for
    -> r -> index(s) for the delete row
    -> c -> index(s) for the delete column

    Note: A MinorMatrix is a smaller matrix obtained from a larget one by deleting one or more rows and/or one columns
    */

    std::vector<Matrix::size_type> mDim = m.dimension();
    std::size_t rSize = r.size();
    std::size_t cSize = c.size();
    Matrix::size_type mmRow = mDim[0] - rSize;
    Matrix::size_type mmCol = mDim[1] - cSize;

    Matrix mm = Matrix(mmRow,mmCol); // Initialize zero matrix we will fill in

    Matrix::size_type curr_mmRow = 0;
    for (Matrix::size_type i = 0; i < mDim[0]; ++i) {
        if (r.find(i) != r.end()) {
            continue;
        }
        Matrix::size_type curr_mmCol = 0;
        for (Matrix::size_type j = 0; j < mDim[1]; ++j) {
            if (c.find(j) != c.end()) {
                continue;
            }
            mm[curr_mmRow][curr_mmCol] = m[i][j];
            ++curr_mmCol;
        }
        ++curr_mmRow;
    }

    return mm;

}

double Determinant(const Matrix& m) {
    std::vector<Matrix::size_type> mDim = m.dimension();
    if (mDim[0] != mDim[1]) {
        throw std::invalid_argument("Matrix must be a square matrix to compute the determinant.");
    }

    if (mDim[0] == 1) {
        return m[0][0]; // 1x1 determinant has matrix equal to its only entry
    }
    else if (mDim[0] == 2) {
        return Determinant2x2(m);
    }
    else { // Use Cofactor Expansion
        Matrix::size_type r_delete = 0;             // row to delete for MinorMatrix
        std::vector<Matrix::size_type> z_indices;   // indicates col index in cofactor row where entry = 0
        Matrix::size_type z_count = 0;                            // indicates # of zeros in cofactor row

        // Below for loop searches for row with the most zeros to use as delete row in MinorMatrix
        for (Matrix::size_type i = 0; i < mDim[0]; ++i) {
            if (z_count == mDim[1]) {
                return 0.0; // If zero matrix, determinant = 0
            }

            Matrix::size_type r_num_zero = 0; // Current # of 0 entries in row
            std::vector<Matrix::size_type> r_z_indices; // col index for current row check 0 entries

            for (Matrix::size_type j = 0; j < mDim[1]; ++j) {
                if (m[i][j] == 0) {
                    ++r_num_zero;
                    r_z_indices.push_back(j);
                }
            }

            if (r_num_zero > z_count) {
                r_delete = i;
                z_count = r_num_zero;
                z_indices = std::move(r_z_indices);
            }
        }

        // Now we search for the first nonzero entry in our row to serve as delete col in MinorMatrix
        double det = 0.0;
        for (Matrix::size_type k = 0; k < mDim[1]; ++k) {
            if (std::find(z_indices.begin(), z_indices.end(), k) != z_indices.end()) {
                continue;
            }
            else {
                Matrix mm = MinorMatrix(m, std::set<Matrix::size_type>{r_delete}, std::set<Matrix::size_type>{k});

                double sign = ((r_delete + k) & 1) ? -1.0 : 1.0;    // bitwise parity
                det += std::pow(-1,(r_delete + k)) * m[r_delete][k] * Determinant(mm);
            }
        }

        return det;

    }
    
}

Matrix RowScale(const Matrix& m, const Matrix::size_type& r) {
    double f = 0.0;
    for (double x : m[r]) {
        if (x != 0) {
            f = x;
            break;
        }
    }

    if (f == 0.0) {
        return m;
    }
    else {
        Matrix m2(m);
    
        for (double& c : m2[r]) {
            c /= f;
        }

        return m2;
    }
}

Matrix RowSwap(const Matrix& m, const Matrix::size_type& r1, const Matrix::size_type& r2) {
    Matrix rs(m);
    std::swap(rs[r1],rs[r2]);
    return rs;
}

Matrix RowReplacement(const Matrix&m, const Matrix::size_type& r1, const Matrix::size_type& r2) {
    /*
    r1 -> Index for target row (one being changed)
    r2 -> Index for source/pivot row (one used to modify target row)
    
    Proper Row Replacement formula is R1 -> R1 + kR2, where k = -x/p -> (p = R2 pivot value, x = R1 value in R2's pivot column)
    */
    
    // If target or source row is zero vector, do nothing
    if (la::isZeroVector(m[r1]) || la::isZeroVector(m[r2])) {
        return m;
    }

    // Find source row's pivot column
    Matrix::size_type sp;
    for (Matrix::size_type i = 0; i < m[r2].size(); ++i) {
        if (m[r2][i] != 0.0) {
            sp = i;
            break;
        }
    }

    // Check target row value at source's pivot column. If 0, we do nothing
    if (m[r1][sp] == 0.0) {
        return m;
    }

    // Calculate k
    double k = -1 * m[r1][sp] / m[r2][sp];

    // Find new target row by 
    std::vector<double> nr1;
    for (Matrix::size_type j = 0; j < m[r1].size(); ++j) {
        double newVal = m[r1][j] + (k * m[r2][j]);
        nr1.push_back(newVal);
    }

    // Initialize return matrix
    la::Matrix rrm(m);
    rrm[r1] = la::Vector(nr1);

    return rrm;
}

Matrix RowEchelonForm(const Matrix& m) {
    /*
    Here are steps for the foward elimination algorithm to find REF:
    ----------------------------------------------------------------
    Step 1: Identify the pivot
         -Start with the first column and top row
        - Find the first nonzero entry in that column
        - Swap rows, if needed, so the first nonzero entry is at the top.
        - If no rows have a nonzero entry in this column, move one column to the right and repeat this step
    Step 2: Create zeros for each element below the pivot
        - Use row replacement to eliminate nonzero entries in the pivot column for all rows below the current pivot row
    Step 3: Move to the next pivot column
        - Move one row down and one column to the right
    Step 4: Repeat
        - Repeat steps 1-3 until you run out of rows or all remaining rows are zeros
    */
    Matrix::size_type rc = m.rowCount();
    Matrix::size_type cc = m.colCount();
    Matrix::size_type curr_r = 0;
    Matrix::size_type curr_c = 0;

    la::Matrix ref(m);
    
    while (curr_r < rc && curr_c < cc) {
        // Check current row/col to see if it's a pivot
        if (ref[curr_r][curr_c] == 0.0) {

            // Iterate through below rows to search for non-zero to swap with
            for (Matrix::size_type i = curr_r + 1; i < rc; ++i) {
                if (ref[i][curr_c] != 0) {
                    ref = RowSwap(ref,curr_r,i);
                    break;
                }
            }

            // If value is still 0, we have a free variable. Move to next column
            if (ref[curr_r][curr_c] == 0.0) {
                curr_c++;
                continue;
            }
        }

        // Assume we have non-zero value, so we  use row replacement to remove entries below the pivot
        for (Matrix::size_type i = curr_r + 1; i < rc; ++i) {
            ref = RowReplacement(ref,i,curr_r);
        }

        // Increment our current row/col
        curr_r++;
        curr_c++;

    }

    return ref;

}

Matrix ReducedRowEchelonForm(const Matrix& m) {
    /*
    Here are steps for the backward elimination algorithm to find RREF:
    -------------------------------------------------------------------
    Step 0:
        - Execute RowEchelonForm() to get matrix in REF
    Step 1:
        - Start at the bottom-most pivot and use row scaling to set the pivot equal to 1 (if it’s not already 1)
    Step 2:
        - Use row replacement to make every entry above the pivot equal to 0
    Step 3:
        - Move to the next pivot above, and repeat the process until all pivots are equal to 1 and isolated in their columns
    */

    Matrix ref = RowEchelonForm(m);
    std::vector<Matrix::size_type> refDim = ref.dimension();
    std::vector<Matrix::size_type> pvts;

    // Find pivots, starting with right-most
    for (int r = refDim[0] - 1; r >= 0; --r) {
        // make r int bc Matrix::size_type is signed, so --0 becomes a long long int. Use static_cast below to access elements
        for (Matrix::size_type c = 0; c < refDim[1]; ++c) {
            if (ref[static_cast<Matrix::size_type>(r)][c] != 0.0) {
                pvts.push_back(static_cast<Matrix::size_type>(r));
                break;
            }
        }
    }

    for (Matrix::size_type p : pvts) {
        ref = RowScale(ref,p); // Scale row to have pivot = 1
        for (Matrix::size_type k = 0; k < p; ++k) {
            ref = RowReplacement(ref,k,p); // replace rows above to have pivot column val = 0
        }
    }

    return ref;

}

Matrix::size_type Rank(const Matrix& m) {
    /*
    Rank of matrix is number of linearly independent rows/columns.
    We can easily find rank by reducing matrix to REF and finding # of pivots
    */

    Matrix ref = RowEchelonForm(m);
    
    Matrix::size_type r = ref.rowCount(); // Start with Rank = rowCount. Decrement below as needed
    for (Vector& v : ref) {
        if (isZeroVector(v)) {
            r--;
        }
    }

    return r;

}

bool isFullRank(const Matrix& m) {
    /*
    Full Rank of a mxn dimension matrix A means -> Rank(A) = min(m,n)
    */
    std::vector<Matrix::size_type> mDim = m.dimension();
    Matrix::size_type minDim = std::min(mDim[0],mDim[1]);
    Matrix::size_type mRank = Rank(m);

    return (minDim == mRank);
}

Matrix AugmentedMatrix(const Matrix& m1, const Matrix& m2) {
    // Need same row count for m1 and m2 to create augmented matrix
    Matrix::size_type m1RowCount = m1.rowCount();
    Matrix::size_type m1ColCount = m1.colCount();
    Matrix::size_type m2RowCount = m2.rowCount();
    Matrix::size_type m2ColCount = m2.colCount();
    if (m1RowCount != m2RowCount) {
        throw std::invalid_argument("Cannot create augmented matrix of matrices with different row dimension.");
    }

    Matrix am;

    for (Matrix::size_type i = 0; i < m1RowCount; ++i) {
        std::vector<double> cv;
        cv.reserve(m1ColCount + m2ColCount);
        const Vector& row1 = m1[i];
        const Vector& row2 = m2[i];
        cv.insert(cv.end(), row1.begin(), row1.end());
        cv.insert(cv.end(), row2.begin(), row2.end());

        Vector amRow(cv);
        am.push_back(amRow);
    }

    return am;

}

Matrix AugmentedMatrixSplit(const Matrix& m, Matrix::size_type r, AmSide side) {
    // Function takes a matrix and splits it at the r-index (exclusive). Return left or right side
    Matrix ams;
    Matrix::size_type b; // begin index (inclusive)
    Matrix::size_type e; // end index (exclusive)

    switch (side) {
        case la::AmSide::L:
            b = 0;
            e = r;
            break;
        case la::AmSide::R:
            b = r;
            e = m.colCount();
            break;
    }
    
    for (Matrix::size_type i = 0; i < m.rowCount(); ++i) {
        std::vector<double> row;
        for (Matrix::size_type j = b; j < e; ++j) {
            row.push_back(m[i][j]);
        }
        ams.push_back(Vector(row));
    }

    return ams;
}

Matrix InverseMatrix(const Matrix & m) {
    /*
    Matrix A has an inverse only if the following are true:
    1. It's a square matrix
    2. Det(A) != 0
    3. Rank(A) is full rank

    If Matrix satisfies above, here's how to find Inverse:
    ------------------------------------------------------
    - 1-dimension matrix
        - A = [a] and A^-1 = [1/a]
    - 2-dimension matrix
        - A = [[a,b],[c,d]] and A^-1 = (1/det(A)) * [[d,-b],[-c,a]]
    - n-dimension matrix (where n > 2)
        - Use Gaussian elimination...
        - [A|I] --RREF--> [I|A^-1]
        - Steps:
            - Create Augmented Matrix with A and I (identity matrix)
            - Reduce Augmented Matrix to RREF
            - Left side becomes I and right side is A^-1
    */
    Matrix::size_type rc = m.rowCount();

    if (rc != m.colCount()) {
        throw std::invalid_argument("Can only find Inverse of a square matrix.");
    }
    else if (Determinant(m) == 0.0) {
        throw std::invalid_argument("Determinant = 0, so matrix is not invertible.");
    }
    else if (!isFullRank(m)) {
        throw std::invalid_argument("Matrix is not full rank, so there is no inverse matrix.");
    }

    Matrix ident(rc); // Initialize Identity matrix
    Matrix am = AugmentedMatrix(m,ident);

    Matrix amRREF = ReducedRowEchelonForm(am);
    
    return AugmentedMatrixSplit(amRREF,rc,la::AmSide::R);
    // return amRREF;
}

Matrix Transpose(const Matrix& m) {
    /*
    Transpose swaps all columns into rows
    */
    Matrix::size_type mRowCount = m.rowCount();
    Matrix::size_type mColCount = m.colCount();
    Matrix t; // Initialize tranpose matrix

    // Loop through each column, create new Vector, add to t
    for (Matrix::size_type i = 0; i < mColCount; ++i) {
        std::vector<double> t_row;
        t_row.reserve(mColCount);

        for (Matrix::size_type j = 0; j < mRowCount; ++j) {
            t_row.push_back(m[j][i]);
        }

        t.push_back(Vector(t_row));
    }

    return t;
}


bool isOrthonormal(const Matrix&m) {
    /*
    Matrix A is orthonormal if:
    - A is square (n-dimension)
    - A^T multiplied by A = I (identity matrix in n-dimension)
    */
   if (!isSquareMatrix(m)) {
       return false;
    }
    
    Matrix mT = Transpose(m);
    Matrix mm = mT * m;
    
    return isIdentityMatrix(mm);
    
}

Vector Diagonal(const Matrix& m) {
    if (!isSquareMatrix(m)) {
        throw std::invalid_argument("Can only find diagonal of a square matrix.");
    }
    Matrix::size_type mRowCount = m.rowCount();
    std::vector<double> d;
    d.reserve(mRowCount);
    
    for (Matrix::size_type i = 0; i < mRowCount; ++i) {
        d.push_back(m[i][i]);
    }
    
    return Vector(d);
}

double Trace(const Matrix& m) {
    Vector d = Diagonal(m);
    double t;
    for (double e : d) {
        t += e;
    }
    
    return t;
}

bool isOrthogonal(const Matrix& m1, const Matrix& m2) {
    /*
    Two Matrices, A and B are orthogonal if their inner product equals 0.
    That is  -> trace(A^T x B) = 0
    */
    Matrix m1Transpose = Transpose(m1);
    Matrix mm = m1Transpose * m2;

    double t = Trace(mm);
    return t == 0.0;
}

bool isSymmetric(const Matrix& m) {
    // Note: we need to add a "==" operator to Matrix class to compare two matrices
    // On above, we also need a "==" operator for the Vector class
    Matrix mT = Transpose(m);

    return m == mT;
}

bool isDiagonal(const Matrix& m) {
    if (!isSquareMatrix(m)) {
        return false;
    }
    std::vector<Matrix::size_type> mDim = m.dimension();
    for (Matrix::size_type i = 0; i < mDim[0]; ++i) {
        for (Matrix::size_type j = 0; j < mDim[1]; ++j) {
            if (i != j && m[i][j] != 0.0) {
                return false;
            }
        }
    }

    return true;
}

bool isTriangular(const Matrix& m) {
    if (!isSquareMatrix(m)) {
        return false;
    }

    Matrix::size_type mRowCount = m.rowCount();
    bool isUpper = true;
    bool isLower = true;

    for (Matrix::size_type i = 0; i < mRowCount; ++i) {
        for (Matrix::size_type j = 0; j < mRowCount; ++j) {
            if (i > j && m[i][j] != 0.0) {
                isUpper = false;
            }
            if (i < j && m[i][j] != 0.0) {
                isLower = false;
            }
        }
    }
    
    return isUpper || isLower;
}

std::vector<std::string> TriangularType(const Matrix& m) {
    if (!isTriangular(m)) {
        throw std::invalid_argument("Matrix is not Triangular.");
    }

    Matrix::size_type mRowCount = m.rowCount();
    bool isUpper = true;
    bool isLower = true;

    std::vector<std::string> tt;

    for (Matrix::size_type i = 0; i < mRowCount; ++i) {
        for (Matrix::size_type j = 0; j < mRowCount; ++j) {
            if (i > j && m[i][j] != 0.0) {
                isUpper = false;
            }
            if (i < j && m[i][j] != 0.0) {
                isLower = false;
            }
        }
    }

    if (isUpper) {
        tt.push_back("Upper");
    }
    
    if (isLower) {
        tt.push_back("Lower");
    }

    return tt;
}

}