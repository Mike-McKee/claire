#pragma once
#include "Vector.h"

namespace la {

    // Arithmetic
    Vector VectorAdddition(const Vector& v1, const Vector& v2);
    Vector VectorSubtraction(const Vector& v1, const Vector& v2);
    Vector ScalarMultiplication(const Vector& v, double k);
    double DotProduct(const Vector& v1, const Vector& v2);
    // Vector CrossProduct(const Vector& v1, const Vector& v2);

    // Magnitude and Direction
    double VectorNorm(const Vector& v);
    double Distance(const Vector& v1, const Vector& v2);
    Vector Normalize(const Vector& v);
    double AngleBetween(const Vector& v1, const Vector& v2);

    // Transformations and Projections
    Vector Projection(const Vector& v1, const Vector& v2);
    // Vector Reflection(const Vector& v1, const Vector& v2);
    // Vector Rotation(const Vector& v, double theta);

    // Helper Functions
    void PrintVector(const Vector& v);
    bool isZeroVector(const Vector& v);
    bool isSameDimension(const Vector& v1, const Vector& v2);

}