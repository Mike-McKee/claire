#include <iostream>
#include "include/claire/la/Vector.h"
#include "include/claire/la/VectorOps.h"
#include "include/claire/la/MatrixOps.h"

// To Compile:

int main() {
    std::cout << std::endl;
    // Need to test all functions in MatrixOps.h

    // -------- MatrixAddition --------
    std::cout << "-------- MatrixAddition --------\n\n";
    la::Vector a1(std::vector<double>{2,-1});
    la::Vector a2(std::vector<double>{4,3});

    la::Vector b1(std::vector<double>{5,2});
    la::Vector b2(std::vector<double>{-3,6});

    la::Matrix a({a1,a2});
    la::Matrix b({b1,b2});

    
    // -------- MatrixSubtraction --------
    std::cout << "-------- MatrixSubtraction --------\n\n";


    // -------- ScalarMultiplication --------
    std::cout << "-------- ScalarMultiplication --------\n\n";


    // -------- MatrixMultiplication --------
    std::cout << "-------- MatrixAddition --------\n\n";

    std::cout << std::endl;
    return 0;
}