#include <iostream>
#include "claire/la/Matrix.h"
#include "claire/la/VectorOps.h"

// To Compile: g++ .\MatrixTest.cpp ..\src\Vector.cpp ..\src\Matrix.cpp ..\src\VectorOps.cpp -I..\include\ -o MatrixTest

int main() {
    // Need to test all Matrix class methods

    // Initialize vectors
    la::Vector u({1,2,3,4,5});
    la::Vector v({1,2,3,4,5});
    la::Vector w({1,2,3,4,5});
    la::Vector x({1,2,3,4,5});

    // -------- Constructors --------

    // Zero Matrix
    la::Matrix zm(static_cast<la::Matrix::size_type>(3),static_cast<la::Matrix::size_type>(3));
    std::cout << "\n------ Zero Matrix Test ------\n\n";
    for (la::Vector& r : zm) {
        la::PrintVector(r);
        std::cout << "\n";
    }


    // Identity Matrix
    la::Matrix im(5);
    std::cout << "\n------ Identity Matrix Test ------\n\n";

    for (la::Vector& r : im) {
        la::PrintVector(r);
        std::cout << "\n";
    }

    // 3D rotation matrix
    la::Matrix rm('x',90.0);
    std::cout << "\n------ 3D Rotation Matrix Test ------\n\n";

    for (la::Vector& r : rm) {
        la::PrintVector(r);
        std::cout << "\n";
    }

    // Matrix with specific values
    la::Matrix m({u,v,w,x});
    std::cout << "\n------ Matrix With Specific Values Test ------\n\n";

    for (la::Vector& r : m) {
        la::PrintVector(r);
        std::cout << "\n";
    }

    // Copy Constructor
    la::Matrix cm(m);
    std::cout << "\n------ Copy Constructor Test ------\n\n";

    for (la::Vector& r : cm) {
        la::PrintVector(r);
        std::cout << "\n";
    }
    
    // -------- Object Info --------
    std::cout << "\n------ Dimension Test ------\n\n";
    std::cout << "For 5x5 Identity matrix above...\n";
    std::cout << im.dimension()[0] << "x" << im.dimension()[1] << "\n";
    
    std::cout << "\n------ rowCount Test ------\n\n";
    std::cout << im.rowCount() << "\n";

    std::cout << "\n------ colCount Test ------\n\n";
    std::cout << im.colCount() << "\n";
    
    // -------- Operators --------
    std::cout << "\nUsing the 3D Rotation Matrix above, let's access and print the second row...\n";
    std::cout << "Read operator -> ";
    const la::Vector& r2 = rm.row(1);
    la::PrintVector(r2);

    std::cout << "\nWrite operator -> ";
    la::Vector& wr2 = rm.row(1);
    la::PrintVector(wr2);

    std::cout << std::endl;
    return 0;
}