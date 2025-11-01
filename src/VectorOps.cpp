#define _USE_MATH_DEFINES
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
// Vector CrossProduct(const Vector& v1, const Vector& v2) {
//     // Need to verify that v1 and v2 have size = 3, since cross product only exists in 3D
//     if (v1.size() != 3) {
//         throw std::invalid_argument("Cannot perform cross product if vectors are not 3-dimensional.");
//     }
//     else if (!isSameDimension(v1,v2)) {
//         throw std::invalid_argument("Cannot perform cross product if vectors are not same size.");
//     }

// }   

double VectorNorm(const Vector& v) {
    double sum = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += v[i] * v[i];
    }
    return std::sqrt(sum);
}

double Distance(const Vector& v1, const Vector& v2) {
    Vector v1Minusv2 = VectorSubtraction(v1,v2);
    return VectorNorm(v1Minusv2);
}

Vector Normalize(const Vector& v) {
    double vNorm = VectorNorm(v);
    if (vNorm == 0) {
        throw std::runtime_error("Division by Zero error encountered");
    }
    double k = 1.0 / vNorm;
    return ScalarMultiplication(v,k);
}

double AngleBetween(const Vector& v1, const Vector& v2) {
    // Angle Theta between vectors v1, v2 is cos(theta) = (v1 dot v2) / (v1 norm * v2 norm)
    // Then use std::acos() to find angle in radians
    // Note: can potentially create bool parameter that allows us to return in radians insted of degrees
    double v1Norm = VectorNorm(v1);
    double v2Norm = VectorNorm(v2);
    double v1Dotv2 = DotProduct(v1,v2);
    double cosTheta = v1Dotv2 / (v1Norm * v2Norm);
    double rad = std::acos(cosTheta);
    
    // Now convert radians to degrees (M_PI represents standard pi)
    return rad * (180.0 / M_PI);
}

Vector Projection(const Vector& v1, const Vector& v2) {
    double v1Dotv2 = DotProduct(v1,v2);
    double v2Dotv2 = DotProduct(v2,v2);
    double k = v1Dotv2 / v2Dotv2;
    return ScalarMultiplication(v2,k);
}

// Vector Reflection(const Vector& v1, const Vector& v2) {

// }

}