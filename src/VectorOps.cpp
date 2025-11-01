#include <stdexcept>
#include <cmath>
#include <iostream>
#include "claire/la/VectorOps.h"

namespace la {

void PrintVector(const Vector& v) {
    const auto n = v.size();
    std::cout << "[";
    for (std::size_t i = 0; i < n; ++i) {
        if (i == n - 1) {
            std::cout << v[i] << "]";
        }
        else {
            std::cout << v[i] << ",";
        }
    }
}

bool isZeroVector(const Vector& v) {
    const auto& s = v.size();
    for (std::size_t i = 0; i < s; ++i) {
        if (v[i] != 0) {
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

Vector VectorAddition(const Vector& v1, const Vector& v2) {
    if (!isSameDimension(v1,v2)) {
        throw std::invalid_argument("Cannot add vectors with different dimensions");
    }
    std::vector<double> va;
    va.reserve(v1.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
        va.push_back(v1[i] + v2[i]);
    }

    return Vector(va);
}

Vector VectorSubtraction(const Vector& v1, const Vector& v2) {
    if (!isSameDimension(v1,v2)) {
        throw std::invalid_argument("Cannot subtract vectors with different dimensions");
    }
    Vector v2Neg = ScalarMultiplication(v2,-1);
    return VectorAddition(v1,v2Neg);
}

Vector ScalarMultiplication(const Vector& v, double k) {
    std::vector<double> sm;
    sm.reserve(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        sm.push_back(v[i] * k);
    }

    return Vector(sm);
}

double DotProduct(const Vector& v1, const Vector& v2) {
    if (!isSameDimension(v1,v2)) {
        throw std::invalid_argument("Cannot find dot product of vectors with different dimensions");
    }
    double dot = 0;
    for (std::size_t i = 0; i < v1.size(); ++i) {
        dot += v1[i] * v2[i];
    }
    return dot;
}

// TODO: BELOW
// Vector CrossProduct(const Vector& v1, const Vector& v2);

}