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

    // -------- isOrthonormal --------
    std::cout << "\n-------- isOrthonormal--------\n\n";
    la::Vector on1_1(std::vector<double>{0,1});
    la::Vector on1_2(std::vector<double>{-1,0});
    la::Matrix on1({on1_1,on1_2});

    std::cout << "Following matrix should be orthonormal...\n";
    la::PrintMatrix(on1);
    std::cout << "Result from isOrthonormal() is " << std::boolalpha << la::isOrthonormal(on1) << "\n";

    // -------- Diagnoal --------
    std::cout << "\n-------- Diagonal--------\n\n";

    std::cout << "Find diagonal of below matrix...\n\n";
    la::PrintMatrix(am1);
    std::cout << "\nThe diagonal of the above is ";
    la::PrintVector(la::Diagonal(am1));
    std::cout << "\n\n";

    // -------- Trace --------
    std::cout << "\n-------- Trace --------\n\n";
    std::cout << "The trace of the above matrix is " << la::Trace(am1) << "\n\n";


    // -------- isOrthogonal --------
    std::cout << "\n-------- isOrthogonal--------\n\n";
    la::Vector og1_1(std::vector<double>{4,3});
    la::Vector og1_2(std::vector<double>{-1,-2});
    la::Vector og2_1(std::vector<double>{5,0});
    la::Vector og2_2(std::vector<double>{0,10});

    la::Matrix og1({og1_1,og1_2});
    la::Matrix og2({og2_1,og2_2});

    std::cout << "isOrthogonal = " << std::boolalpha << la::isOrthogonal(og1,og2) << "\n";

    // -------- isSymmetric --------
    std::cout << "\n-------- isSymmetric --------\n\n";
    la::Vector sym1(std::vector<double>{1,2,3});
    la::Vector sym2(std::vector<double>{2,0,4});
    la::Vector sym3(std::vector<double>{3,4,1});
    la::Matrix sym({sym1,sym2,sym3});
    std::cout << "isSymmetric = " << std::boolalpha << la::isSymmetric(sym) << "\n";

    // -------- isTriangular --------
    std::cout << "\n-------- isTriangular --------\n\n";
    la::Vector tri1(std::vector<double>{1,0,0,0});
    la::Vector tri2(std::vector<double>{0,2,0,0});
    la::Vector tri3(std::vector<double>{0,0,3,0});
    la::Vector tri4(std::vector<double>{0,0,0,4});
    la::Matrix tri({tri1,tri2,tri3,tri4});
    std::cout << "isTriangular = " << std::boolalpha << la::isTriangular(tri) << "\n";

    // -------- TriangularType --------
    std::cout << "\n-------- TriangularType --------\n\n";
    std::vector<std::string> tt = la::TriangularType(tri);
    for (std::string typ : tt) {
        std::cout << typ << ", ";
    }

    // -------- isDiagonal --------
    std::cout << "\n-------- isDiagonal --------\n\n";
    la::Vector diag1(std::vector<double>{0,1,0,0});
    la::Vector diag2(std::vector<double>{0,0,0,0});
    la::Vector diag3(std::vector<double>{0,0,0,0});
    la::Vector diag4(std::vector<double>{0,0,0,0});
    la::Matrix diag({diag1,diag2,diag3,diag4});
    std::cout << "isDiagonal = " << std::boolalpha << la::isDiagonal(diag) << "\n";

    std::cout << std::endl;
    return 0;
}
