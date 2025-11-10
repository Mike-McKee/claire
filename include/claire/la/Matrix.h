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
    // bool isZeroVector() const;

    //-------- Operators --------
    const Vector& row(size_type i) const;
    Vector& row(size_type i);
    Matrix& operator=(const Matrix& other) {
        if (this == &other) {return *this;}

        rows_ = other.rows_;
        return *this;
    }
};

}