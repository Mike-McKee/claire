#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/MatrixOps.h"

// To Compile: g++ .\MatrixOpsTest3.cpp ..\src\Vector.cpp ..\src\VectorOps.cpp ..\src\Matrix.cpp ..\src\MatrixOps.cpp -I..\include\ -o MatrixOpsTest3

int main() {
    std::cout << std::endl;
    /*
    Functions to test here:
    -----------------------
    SUCCESS:
    Matrix::size_type Rank(const Matrix& m);
    bool isFullRank(const Matrix& m);
    Matrix AugmentedMatrix(const Matrix& m1, const Matrix& m2);
    Matrix AugmentedMatrixSplit(const Matrix& m, Matrix::size_type r, AmSide side);
    Matrix InverseMatrix(const Matrix & m);
    Matrix Transpose(const Matrix& m);
    
    TODO:
    bool isOrthogonal(const Matrix& m1, const Matrix& m2);
    bool isOrthonormal(const Matrix&m);
    Vector Diagonal(const Matrix& m);
    double Trace(const Matrix& m);
    bool isSymmetric(const Matrix& m);
    bool isDiagonal(const Matrix& m);
    bool isTriangular(const Matrix& m);
    std::vector<TriangularMatrixType> TriangularType(const Matrix& m);
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

    // -------- AugmentedMatrix and AugmentedMatrixSplit --------
    std::cout << "-------- AugmentedMatrix and AugmentedMatrixSplit --------\n\n";
    la::Vector am1_1(std::vector<double>({1,2,3}));
    la::Vector am1_2(std::vector<double>({4,5,6}));
    la::Vector am1_3(std::vector<double>({7,8,9}));

    la::Vector am2_1(std::vector<double>({1,0,0}));
    la::Vector am2_2(std::vector<double>({0,1,0}));
    la::Vector am2_3(std::vector<double>({0,0,1}));

    la::Matrix am1({am1_1,am1_2,am1_3});
    la::Matrix am2({am2_1,am2_2,am2_3});

    la::Matrix am1_am2 = la::AugmentedMatrix(am1,am2);
    la::PrintMatrix(am1_am2);
    la::Matrix amsplit = la::AugmentedMatrixSplit(am1_am2,3,la::AmSide::R);
    la::PrintMatrix(amsplit);

    // -------- InverseMatrix --------
    std::cout << "\n-------- InverseMatrix --------\n\n";
    la::Vector im1(std::vector<double>({2,1}));
    la::Vector im2(std::vector<double>({5,3}));
    la::Matrix im({im1,im2});
    la::PrintMatrix(im);
    std::cout << "\nInverse of above matrix is...\n";
    la::Matrix inverse = la::InverseMatrix(im);
    la::PrintMatrix(inverse);

    // -------- Transpose --------
    std::cout << "\n-------- Transpose --------\n\n";
    std::cout << "Transpose of\n";
    la::PrintMatrix(am1);
    std::cout << "is\n";
    la::Matrix am1Transpose = la::Transpose(am1);
    la::PrintMatrix(am1Transpose);

    // -------- isOrthogonal --------
    std::cout << "\n-------- isOrthogonal--------\n\n";
    
    std::cout << std::endl;
    return 0;
}
