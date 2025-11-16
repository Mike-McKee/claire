#define _USE_MATH_DEFINES
// #include <numbers>
#include <stdexcept>
#include <utility>
#include <cmath>
#include "claire/la/VectorOps.h"
#include "claire/la/Matrix.h"

namespace la {

Matrix::Matrix(std::initializer_list<Vector> rows) : rows_(std::move(rows)) {}

Matrix::Matrix(Matrix::size_type rows, Matrix::size_type cols) {
    if (rows == 0 || cols == 0) {
        throw std::invalid_argument("Matrix dimensions must be > 0.");
    }
    rows_.reserve(rows);

    for (Matrix::size_type i = 0; i < rows; ++i) {
        rows_.emplace_back(Vector(cols,0.0));
    }
}

Matrix::Matrix(Matrix::size_type n) {
    if (n == 0) {
        throw std::invalid_argument("Matrix dimensions must be > 0.");
    }

    std::vector<Vector> rows;
    rows.reserve(n);

    for (Matrix::size_type i = 0; i < n; ++i) {
        Vector r(n, 0.0);
        r[i] = 1.0;
        rows.push_back(r);
    }

    rows_ = std::move(rows);
}

Matrix::Matrix(char axis, double degrees) {
    // Note: This only works for 3D matrix
    // Matrix is used as transformation matrix over x, y, or z axis

    // Verify axis = 'x', 'y', or 'z'
    if (axis != 'x' && axis != 'y' && axis != 'z') {
        throw std::invalid_argument("Rotation axis must be x, y, or z.");
    }
    
    Vector r1 = Vector(3,0.0);
    Vector r2 = Vector(3,0.0);
    Vector r3 = Vector(3,0.0);

    double c = std::cos(degrees * (M_PI / 180.0));
    double s = std::sin(degrees * (M_PI / 180.0));

    if (axis == 'x') {
        r1[0] = 1.0;
        r2[1] = c;
        r2[2] = -s;
        r3[1] = s;
        r3[2] = c;
    }
    else if (axis == 'y') {
        r1[0] = c;
        r1[2] = s;
        r2[1] = 1.0;
        r3[0] = -s;
        r3[2] = c;
    }
    else if (axis == 'z') {
        r1[0] = c;
        r1[1] = -s;
        r2[0] = s;
        r2[1] = c;
        r3[2] = 1.0;
    }

    rows_ = {std::move(r1), std::move(r2), std::move(r3)};

}

Matrix::Matrix(const std::vector<Vector>& rows) : rows_(rows) {
    if (rows_.empty()) {
        throw std::invalid_argument("Matrix cannot be built with no rows.");
    }

    // Verify each row has the same dimension
    Vector::size_type colCount = rows_[0].size();
    for (const auto& row : rows_) {
        if (row.size() != colCount) {
            throw std::invalid_argument("All rows must have same dimension.");
        }
    }
}

// Matrix::Matrix(const Matrix& other) = default;

const std::vector<Vector>& Matrix::data() const {
    return rows_;
}

Matrix::size_type Matrix::rowCount() const {
    return rows_.size();
}

Matrix::size_type Matrix::colCount() const {
    return rows_[0].size();
}

std::vector<Matrix::size_type> Matrix::dimension() const {
    return {Matrix::rowCount(),Matrix::colCount()};
}

void Matrix::push_back(const Vector& v) {
    if (rows_.empty()) {
        rows_.push_back(v);
    }
    else {
        // verify v has same number of columns as other rows in matrix
        if (v.size() != rows_[0].size()) {
            throw std::invalid_argument("All rows in a matrix must have the same number of columns.");
        }
        rows_.push_back(v);
    }
}

void Matrix::reserve(Matrix::size_type s) {
    rows_.reserve(s);
}

}