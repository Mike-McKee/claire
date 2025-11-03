#define _USE_MATH_DEFINES
#include <stdexcept>
#include <cmath>
#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/MatrixOps.h"

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
Vector CrossProduct(const Vector& v1, const Vector& v2) {
    // Need to verify that v1 and v2 have size = 3, since cross product only exists in 3D
    if (v1.size() != 3) {
        throw std::invalid_argument("Cannot perform cross product if vectors are not 3-dimensional.");
    }
    else if (!isSameDimension(v1,v2)) {
        throw std::invalid_argument("Cannot perform cross product if vectors are not same size.");
    }
    
    double i = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    double j = (v1[0] * v2[2]) - (v1[2] * v2[0]);
    double k = (v1[0] * v2[1]) - (v1[1] * v2[0]);

    return Vector({i,j,k});

}   

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

Vector Rotation2D(const Vector& v, double degrees) {
    // Typically use transformation matrix, but for now I'll use the basic formula:
    // For vector v = [x,y], rotation by theta degress results in rotation vector Rv, where:
    // Rv = [x*cos(theta) - y*sin(theta), x*sin(theta) + y*cos(theta)]
    if (v.size() != 2) {
        throw std::invalid_argument("Can only use Rotation2D() for two-dimensional vectors.");
    }
    
    // double rad = degrees * (M_PI / 180.0);
    // double new_x = (v[0] * cos(rad)) - (v[1] * sin(rad));
    // double new_y = (v[0] * sin(rad)) + (v[1] * cos(rad));

    double rad = degrees * (M_PI / 180.0);
    Vector xrv = Vector(std::vector<double>{cos(rad),-sin(rad)});
    Vector yrv = Vector(std::vector<double>{sin(rad),cos(rad)});

    Matrix rotationMatrix = Matrix(std::vector<Vector>{xrv,yrv});
    Matrix vm = Matrix({v});

    Matrix rotatedv = MatrixMultiplication(rotationMatrix,vm);

    return rotatedv.row(0);
}

Vector Rotation3D(const Vector& v, double degrees, char axis) {
    if (v.size() != 3) {
        throw std::invalid_argument("Must enter a vector in 3-dimension to rotate.");
    }

    if (axis != 'x' && axis != 'y' && axis != 'z') {
        throw std::invalid_argument("Rotation axis must be x, y, or z.");
    }

    Matrix vm = Matrix({v});
    Matrix rotationMatrix = Matrix(3.0,3.0);

    switch (axis) {
        case 'x': {
            rotationMatrix = Matrix('x',degrees);
            break;
        }
        case 'y': {
            rotationMatrix = Matrix('y',degrees);
            break;
        }
        case 'z': {
            rotationMatrix = Matrix('z',degrees);
            break;
        }
    }

    Matrix rotatedVector = MatrixMultiplication(rotationMatrix,vm);

    return rotatedVector.row(0);

}

Vector Reflection2D(const Vector& v, char rtype, double lineslope = 0.0,
                    double xcoeff = 0.0, double ycoeff = 0.0, double intercept = 0.0) {
    /*
    Types of 2D reflections include:
        - about x-axis
        - about y-axis
        - along origin
        - about x = y line

    char rtype -> values = ('x','y','o','l')
    */

    if (rtype != 'x' && rtype != 'y' && rtype != 'o' && rtype != 'l') {
        throw std::invalid_argument("Reflection type must be x, y, o, or l.");
    }

    Vector xReflect = Vector(std::vector<double>{1.0,0.0});
    Vector yReflect = Vector(std::vector<double>{0.0,1.0});
    Matrix reflectionM = Matrix(std::vector<Vector>{xReflect,yReflect});

    switch (rtype) {
        case 'x': {
            reflectionM.row(1)[1] = -1.0;
            break;
        }
        case 'y': {
            reflectionM.row(0)[0] = -1.0;
            break;
        }
        case 'o': {
            reflectionM.row(0)[0] = -1.0;
            reflectionM.row(1)[1] = -1.0;
            break;
        }
        case 'l': {
            /*
            Using Householder transformation here: v' = (I - 2nn^T)x + 2dn
            Where:
            - v = original vector
            - v' = reflected vector
            - I = identity matrix in 2D
            - n = unit normal vector to the line (length = 1 and perpendicular to the line)
            - n^T = transpose of n
            - nn^T = "outer product" = projection matrix of any vector onto normal direction
                - determines how much of x lies along the normal direction
            - d = signed distance from origin to the line
            - I - 2nn^T = Householder reflection matrix = 
            - 2dn = translation correction that adjusts for lines that do not pass through origin

            Input requires line in form -> Ax + By + C = 0
            - xcoeff = A
            - ycoeff = B
            - intercept = C

            - Using A, B, C, we can derive the following equation for reflection:
                v' = (x',y')
                x' = ((BB - AA) / (AA + BB))x - ((2AB)/(AA + BB))y - ((2AC)/(AA + BB))
                y' = -((2AB)/ (AA + BB))x + ((AA - BB)/(AA + BB))y - ((2BC)/(AA + BB))
            */
           
            double a2 = xcoeff * xcoeff;
            double b2 = ycoeff * ycoeff;
            double ab = xcoeff * ycoeff;
            double ac = xcoeff * intercept;
            double bc = ycoeff * intercept;
            
            double x_ref = (v[0] * (b2-a2)/(a2+b2)) - (v[1] * (2*ab)/(a2+b2)) - ((2*ac)/(a2+b2));
            double y_ref = (-1 * v[0] * (2*ab)/(a2+b2)) + (v[1] * (a2-b2)/(a2+b2)) - ((2*bc)/(a2+b2));

            Vector reflectV = Vector(std::vector<double>{x_ref,y_ref});
            
            return reflectV;
        }
    }

    Matrix vm = Matrix({v});
    Matrix reflectedV = MatrixMultiplication(reflectionM,vm);

    return vm.row(0);
}
}