#pragma once
#include <vector>
#include <cstddef>
#include <tuple>
#include "Vector.h"

namespace la {

class Matrix
{
private:
    std::vector<Vector> rows_;
    
public:
    using size_type = std::size_t;

    // -------- Constructor --------
    Matrix() = default; // Empty matrix
    Matrix(std::initializer_list<Vector> rows);
    explicit Matrix(size_type rows, size_type cols); // Zero Matrix
    Matrix(size_type n); // Identity Matrix
    explicit Matrix(char axis, double degrees); // Create a 3D rotation matrix
    Matrix(const std::vector<Vector>& rows); // Matrix with specific values
    Matrix(const Matrix& other) = default; // Copy constructor

    // ~Matrix();

    
    // -------- Accessors --------
    const std::vector<Vector>& data() const;
    std::vector<Vector>::iterator begin() {return rows_.begin();}
    std::vector<Vector>::iterator end() {return rows_.end();}
    std::vector<Vector>::const_iterator begin() const {return rows_.begin();}
    std::vector<Vector>::const_iterator end() const {return rows_.end();}
    // std::vector<Vector>::reverse_iterator rbegin() const {return rows_.rbegin();}
    // std::vector<Vector>::reverse_iterator rend() const {return rows_.rend();}
    
    // -------- Object Info --------
    std::vector<size_type> dimension() const;
    size_type rowCount() const;
    size_type colCount() const;
    
    void push_back(const Vector& v);
    void reserve(size_type s);

    //-------- Operators --------
    Vector operator[](size_type i) const { return rows_[i]; }
    Vector& operator[](size_type i) { return rows_[i]; }
    Matrix& operator=(const Matrix& other) {
        if (this == &other) {return *this;}

        rows_ = other.rows_;
        return *this;
    }
    bool operator==(const Matrix& other) const {
        return this->rows_ == other.rows_;
    }
    Matrix& operator+=(const Matrix& other) {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Matrices must have same dimension for Matrix Addition.");
        }
        for (size_type i = 0; i < rowCount(); ++i) {
            for (size_type j = 0; j < colCount(); ++j) {
                rows_[i][j] += other[i][j];
            }
        }

        return *this;
    }
    Matrix& operator-=(const Matrix& other) {
        if (dimension() != other.dimension()) {
            throw std::invalid_argument("Matrices must have same dimension for Matrix Addition.");
        }
        for (size_type i = 0; i < rowCount(); ++i) {
            for (size_type j = 0; j < colCount(); ++j) {
                rows_[i][j] -= other[i][j];
            }
        }

        return *this;
    }
    Matrix& operator*=(const double& k) {
        for (size_type i = 0; i < rowCount(); ++i) {
            for (size_type j = 0; j < colCount(); ++j) {
                rows_[i][j] *= k;
            }
        }

        return *this;
    }
    Matrix& operator*=(const Matrix& other) {
        if (colCount() != other.rowCount()) {
            throw std::invalid_argument("Left matrix column count must equal right matrix row count for Matrix Multiplication.");
        }
        std::vector<size_type> this_dim = dimension();
        std::vector<size_type> other_dim = other.dimension();
        
        std::vector<Vector> new_this;
        new_this.reserve(this_dim[0]);

        for (size_type i = 0; i < this_dim[0]; ++i) {
            Vector r;
            r.reserve(other_dim[1]);
            for (size_type j = 0; j < other_dim[1]; ++j) {
                double entry = 0.0;
                for (size_type k = 0; k < this_dim[1]; ++k) {
                    entry += rows_[i][k] * other[k][j];
                }
                r.push_back(entry); //need to add to vector class
            }
            new_this.push_back(r);
        }

        rows_ = std::move(new_this);
        return *this;
    }
};

Matrix operator+(Matrix a, const Matrix& b) {
    if (a.dimension() != b.dimension()) {
        throw std::invalid_argument("Matrices must have same dimension for Matrix Addition.");
    }
    for (Matrix::size_type i = 0; i < a.rowCount(); ++i) {
        for (Matrix::size_type j = 0; j < a.colCount(); ++j) {
            a[i][j] += b[i][j];
        }
    }
    return a;
}

Matrix operator-(Matrix a, const Matrix& b) {
    if (a.dimension() != b.dimension()) {
        throw std::invalid_argument("Matrices must have same dimension for Matrix Subtraction.");
    }
    for (Matrix::size_type i = 0; i < a.rowCount(); ++i) {
        for (Matrix::size_type j = 0; j < a.colCount(); ++j) {
            a[i][j] -= b[i][j];
        }
    }
    return a;
}

Matrix operator*(Matrix a, const double& k) {
    for (Matrix::size_type i = 0; i < a.rowCount(); ++i) {
        for (Matrix::size_type j = 0; j < a.colCount(); ++j) {
            a[i][j] *= k;
        }
    }

    return a;
}

Matrix operator*(Matrix a, const Matrix& b) {
    if (a.colCount() != b.rowCount()) {
        throw std::invalid_argument("Left matrix column count must equal right matrix row count for Matrix Multiplication.");
    }
    std::vector<Matrix::size_type> a_dim = a.dimension();
    std::vector<Matrix::size_type> b_dim = b.dimension();

    Matrix m;
    m.reserve(a_dim[0]);

    for (Matrix::size_type i = 0; i < a_dim[0]; ++i) {
        Vector r;
        r.reserve(b_dim[1]);
        for (Matrix::size_type j = 0; j < b_dim[1]; ++j) {
            double entry = 0.0;
            for (Matrix::size_type k = 0; k < a_dim[1]; ++k) {
                entry += a[i][k] *= b[k][j];
            }
            r.push_back(entry);
        }
        m.push_back(r);
    }

    a = m;
    return a;
}

}
