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
    // Matrix(); // Empty matrix

    explicit Matrix(size_type rows, size_type cols); // Zero Matrix
    Matrix(size_type n); // Identity Matrix
    Matrix(char axis, double degrees); // Create a 3D rotation matrix
    Matrix(const std::vector<Vector>& rows); // Matrix with specific values
    Matrix(const Matrix& other) = default; // Copy constructor


    // ~Matrix();

    
    // -------- Accessors --------
    const std::vector<Vector>& data() const;
    
    // -------- Object Info --------
    std::vector<size_type> dimension() const;
    size_type rowCount() const;
    size_type colCount() const;
    // bool isZeroVector() const;

    //-------- Operators --------
    const Vector& row(size_type i) const;
    Vector& row(size_type i);
}

}