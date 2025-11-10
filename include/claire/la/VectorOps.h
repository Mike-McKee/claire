#pragma once
#include "Vector.h"
/*
FUTURE TODO:
------------
Consider making functions directly modify the input vector.
Current method of creating multiple copy-vector objects within functions can become expensive.
*/

namespace la {

    // Arithmetic
    Vector VectorAddition(const Vector& v1, const Vector& v2);
    Vector VectorSubtraction(const Vector& v1, const Vector& v2);
    Vector ScalarMultiplication(const Vector& v, double k);
    double DotProduct(const Vector& v1, const Vector& v2);
    Vector CrossProduct(const Vector& v1, const Vector& v2);

    // Magnitude and Direction
    double VectorNorm(const Vector& v);
    double Distance(const Vector& v1, const Vector& v2);
    Vector Normalize(const Vector& v);
    double AngleBetween(const Vector& v1, const Vector& v2);

    // Transformations and Projections
    Vector Projection(const Vector& v1, const Vector& v2);
    // Vector Rotation2D(const Vector& v, double degrees = 90.0);
    // Vector Rotation3D(const Vector& v, double degrees, char axis);
    // Vector Reflection2D(const Vector& v, char rtype, double lineslope = 0.0,
    //                     double xcoeff = 0.0, double ycoeff = 0.0, double intercept = 0.0);

    // Helper Functions
    void PrintVector(const Vector& v);
    bool isZeroVector(const Vector& v);
    bool isSameDimension(const Vector& v1, const Vector& v2);

}