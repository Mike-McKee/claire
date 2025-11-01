#include <stdexcept>
#include <cmath>
#include <iostream>
#include "claire/la/VectorOps.h"

namespace la {

void PrintVector(const Vector& v) {
    const auto n = v.size();
    const auto& d = v.data();
    std::cout << "[";
    for (std::size_t i = 0; i < n; ++i) {
        if (i == n - 1) {
            std::cout << d[i] << "]";
        }
        else {
            std::cout << d[i] << ",";
        }
    }
}

bool isZeroVector(const Vector& v) {
    const auto& d = v.data();
    const auto& s = v.size();
    for (std::size_t i = 0; i < s; ++i) {
        if (d[i] != 0) {
            return false;
        }
    }
    return true;
}

bool isSameDimension(const Vector& v1, const Vector& v2) {
    if (v1.size() == v2.size()) {
        return true;
    }
    else {
        return false;
    }
}

Vector VectorAdddition(const Vector& v1, const Vector& v2) {
    if (!isSameDimension(v1,v2)) {
        throw std::invalid_argument("Cannot add vectors with different dimensions");
    }
    const auto& a = v1.data();
    const auto& b = v2.data();
    std::vector<double> va;
    va.reserve(5);

    for (std::size_t i = 0; i < a.size(); ++i) {
        va.push_back(a[i] + b[i]);
    }

    return Vector(va);
}

}