#include <iostream>
// #include "include/claire/la/Vector.h"
#include "claire/la/VectorOps.h"
#include "claire/la/MatrixOps.h"

// To Compile: g++ .\MatrixOpsTest.cpp ..\src\Vector.cpp ..\src\VectorOps.cpp ..\src\Matrix.cpp ..\src\MatrixOps.cpp -I..\include\ -o MatrixOpsTest

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

    std::cout << "Matrix a = \n";
    la::PrintMatrix(a);
    std::cout << "\nMatrix b = \n";
    la::PrintMatrix(b);

    la::Matrix aPlusb = la::MatrixAddition(a,b);

    std::cout << "\na + b = \n";
    la::PrintMatrix(aPlusb);

    // -------- MatrixSubtraction --------
    std::cout << "\n-------- MatrixSubtraction --------\n\n";
    la::Vector c1(std::vector<double>{6,3});
    la::Vector c2(std::vector<double>{2,5});

    la::Vector d1(std::vector<double>{1,4});
    la::Vector d2(std::vector<double>{7,0});

    la::Matrix c({c1,c2});
    la::Matrix d({d1,d2});

    std::cout << "Matrix c = \n";
    la::PrintMatrix(c);
    std::cout << "\nMatrix d = \n";
    la::PrintMatrix(d);

    la::Matrix cMinusd = la::MatrixSubtraction(c,d);

    std::cout << "\nc - d = \n";
    la::PrintMatrix(cMinusd);

    // -------- ScalarMultiplication --------
    std::cout << "\n-------- ScalarMultiplication --------\n\n";
    
    la::Matrix cscalar3 = la::ScalarMultiplication(c,3);
    std::cout << "Matrix c = \n";
    la::PrintMatrix(c);

    std::cout << "\nc * 3 = \n";
    la::PrintMatrix(cscalar3);

    // -------- MatrixMultiplication --------
    std::cout << "\n-------- MatrixMultiplication --------\n\n";
    la::Vector m1(std::vector<double>{2,3});
    la::Vector m2(std::vector<double>{1,4});
    
    la::Vector n1(std::vector<double>{5,2});
    la::Vector n2(std::vector<double>{7,1});

    la::Matrix m({m1,m2});
    la::Matrix n({n1,n2});

    std::cout << "Matrix m = \n";
    la::PrintMatrix(m);
    std::cout << "\nMatrix n = \n";
    la::PrintMatrix(n);

    la::Matrix mtimen = la::MatrixMultiplication(m,n);
    std::cout << "\nm x n = \n";
    la::PrintMatrix(mtimen);
    
    // -------- Determinant2x2 --------
    std::cout << "\n-------- Determinant2x2 --------\n\n";

    double det = la::Determinant2x2(m);

    std::cout << "Expect det(m) = 5\n";
    std::cout << "det(m) = " << det << "\n";


    // -------- Determinant --------
    std::cout << "\n-------- Determinant --------\n\n";

    la::Vector det1(std::vector<double>{2,-1,3});
    la::Vector det2(std::vector<double>{0,4,5});
    la::Vector det3(std::vector<double>{1,-2,1});
    la::Matrix det_1({det1,det2,det3});
    std::cout << "Determinant of matrix -> \n";
    la::PrintMatrix(det_1);
    std::cout << "is " << la::Determinant(det_1) << "\n";

    la::Vector det4(std::vector<double>{1,2,0,-1});
    la::Vector det5(std::vector<double>{3,0,1,2});
    la::Vector det6(std::vector<double>{2,-1,3,0});
    la::Vector det7(std::vector<double>{1,4,-2,1});
    la::Matrix det_2({det4,det5,det6,det7});
    std::cout << "\nDeterminant of matrix -> \n";
    la::PrintMatrix(det_2);
    std::cout << "is " << la::Determinant(det_2) << "\n";

    std::cout << std::endl;
    return 0;
}