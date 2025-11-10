#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/MatrixOps.h"

// To Compile: g++ .\MatrixOpsTest3.cpp ..\src\Vector.cpp ..\src\VectorOps.cpp ..\src\Matrix.cpp ..\src\MatrixOps.cpp -I..\include\ -o MatrixOpsTest3

int main() {
    std::cout << std::endl;
    /*
    Functions to test here:
    -----------------------
    Matrix::size_type Rank(const Matrix& m);
    bool isFullRank(const Matrix& m);
    Matrix InverseMatrix(const Matrix & m);
    Matrix Transpose(const Matrix& m);
    bool isOrthogonal(const Matrix& m1, const Matrix& m2);
    Vector Diagonal(const Matrix& m);
    double Trace(const Matrix& m);
    bool isSymmetric(const Matrix& m);
    bool isDiagonal(const Matrix& m);
    bool isTriangular(const Matrix& m);
    std::string TriangularType(const Matrix& m);
    */

    // -------- Rank and isFullRank --------
    std::cout << "-------- Rank and isFullRank --------\n\n";
    
    la::Vector r1(std::vector<double>({0,2,-3,1}));
    la::Vector r2(std::vector<double>({4,-2,1,7}));
    la::Vector r3(std::vector<double>({-2,1,5,-8}));
    la::Matrix r({r1,r2,r3});
    la::PrintMatrix(r);
    std::cout << "\n\nRank of above matrix = " << la::Rank(r) << "\n";
    std::cout << "Is Full Rank = " << std::boolalpha << la::isFullRank(r) << "\n\n";

    la::Vector s1(std::vector<double>({1,2,3}));
    la::Vector s2(std::vector<double>({2,4,6}));
    la::Vector s3(std::vector<double>({1,1,1}));
    la::Matrix s({s1,s2,s3});
    la::PrintMatrix(s);
    std::cout << "\nRank of above matrix = " << la::Rank(s) << "\n";
    std::cout << "Is Full Rank = " << std::boolalpha << la::isFullRank(s) << "\n";

    std::cout << std::endl;
    return 0;
}
