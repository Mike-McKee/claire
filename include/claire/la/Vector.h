#pragma once
#include <stdexcept>
#include <vector>
#include <cstddef> // size_t

/*
TODO:
- Add a reserve() function that allows us to reserve size for vector
    - This will help for certain operations, especially regarding Vectors in matrices
- Add a push_back() function to add elements to a vector (similar to std::vector push_back())
*/

namespace la {

class Vector
{
private:
    std::vector<double> data_;
public:
    // using value_type = double; May use for future
    using size_type = std::size_t;

    //-------- Constructor --------
    Vector() = default;
    Vector(std::initializer_list<double> data);
    explicit Vector(std::vector<double> v);
    Vector(size_type n, double value); // Creates n-dimension vector with same value for each entry
    // ~Vector();

    //-------- Accessors --------
    const std::vector<double>& data() const;
    std::vector<double>::iterator begin() {return data_.begin();}
    std::vector<double>::iterator end() {return data_.end();}
    std::vector<double>::const_iterator begin() const { return data_.begin(); }
    std::vector<double>::const_iterator end()   const { return data_.end(); }
    // std::vector<double>::reverse_iterator rbegin() {return data_.rbegin();}
    // std::vector<double>::reverse_iterator rend() {return data_.rend();}

    //-------- Object Info --------
    size_type size() const;
    bool isEmpty() const;

    void push_back(const double& d);
    void reserve(size_type s);

    //-------- Operators --------
    Vector& operator+=(const Vector& other) {
        if (size() != other.size()) {
            throw std::invalid_argument("Vectors must be same size for vector addition.");
        }
        for (size_type i = 0; i < size(); ++i) {
            data_[i] += other[i];
        }
        return *this;
    }
    Vector& operator-=(const Vector& other) {
        if (size() != other.size()) {
            throw std::invalid_argument("Vectors must be the same size for vector subtraction.");
        }
        for (size_type i = 0; i < size(); ++i) {
            data_[i] -= other[i];
        }
        return *this;
    }
    Vector& operator*=(const double& k) {
        for (size_type i =0; i < size(); ++i) {
            data_[i] *= k;
        }
        return *this;
    }
    Vector& operator=(const Vector& other) {
        data_ = std::move(other.data_);
        return *this;
    }
    double operator[](size_type i) const { return data_[i]; }
    double& operator[](size_type i) { return data_[i]; }
    bool operator==(const Vector& other) const {
        return this->data_ == other.data_;
    }
};

Vector operator+(Vector u, const Vector& v) {
    if (u.size() != v.size()) {
        throw std::invalid_argument("Vectors must be the same size for vector addition.");
    }
    for (Vector::size_type i = 0; i < u.size(); ++i) {
        u[i] += v[i];
    }

    return u;
}

Vector operator-(Vector u, const Vector& v) {
    if (u.size() != v.size()) {
        throw std::invalid_argument("Vectors must be the same size for vector addition.");
    }
    for (Vector::size_type i = 0; i < u.size(); ++i) {
        u[i] -= v[i];
    }

    return u;
}

Vector operator*(Vector u, const double& k) {
    for (Vector::size_type i = 0; i < u.size(); ++i) {
        u[i] *= k;
    }

    return u;
}

}
