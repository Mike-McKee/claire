#include <iostream>
#include "claire/la/Vector.h"
#include "claire/la/VectorOps.h"

// To Compile: g++ .\VectorTest.cpp ..\src\Vector.cpp ..\src\VectorOps.cpp ..\src\Matrix.cpp -I..\include\  -std=c++20 -o VectorTest

int main() {

    // Here we'll test all Vector class methods

    // Explicit vector constructor
    la::Vector v({1.0,2.0,4.0,5.0,6.3});

    // Create n-dimension vector with same value for each entry
    la::Vector w(5,2.0);

    // Explicit vector constructor with different initialization method
    la::Vector u = la::Vector(5,4.0);

    std::cout << "v = ";
    la::PrintVector(v);
    std::cout << "\n";

    std::cout << "w = ";
    la::PrintVector(w);
    std::cout << "\n";

    std::cout << "u = ";
    la::PrintVector(u);
    std::cout << "\n";
    
    // Test size()
    std::cout << "v size = " << v.size() << "\n";
    std::cout << "w size = " << w.size() << "\n";
    std::cout << "u size = " << u.size() << "\n";

    la::Vector x(0,0.0);
    std::cout << "is v empty = " << v.isEmpty() << "\n";
    std::cout << "is x empty = " << x.isEmpty() << "\n";

    std::cout << "v[0] = " << v[0] << "\n";
    std::cout << "w[3] = " << w[0] << "\n";
    std::cout << "u[2] = " << u[0] << "\n";
    return 0;
}