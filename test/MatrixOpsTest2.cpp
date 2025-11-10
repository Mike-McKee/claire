#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/MatrixOps.h"

// To Compile: g++ .\MatrixOpsTest2.cpp ..\src\Vector.cpp ..\src\VectorOps.cpp ..\src\Matrix.cpp ..\src\MatrixOps.cpp -I..\include\ -o MatrixOpsTest2

// This script is to test our REF and RREF functions (as well as helper functions)

int main() {
    std::cout << std::endl;

    // -------- RowScale --------
    std::cout << "-------- RowScale --------\n\n";

    la::Vector r1(std::vector<double>{2,-1,3,5});
    la::Vector r2(std::vector<double>{4,0,6,8});
    la::Vector r3(std::vector<double>{-2,5,-9,-3});
    la::Matrix r({r1,r2,r3});
    std::cout << "r is -> \n";
    la::PrintMatrix(r);

    std::cout << "\nLet's scale row 1...\n";
    la::Matrix r1Scaled = la::RowScale(r,0);
    la::PrintMatrix(r1Scaled);
    
    // -------- RowScale --------
    std::cout << "-------- RowScale --------\n\n";
    std::cout << "Swap rows 1 and 2 of r ->\n";
    la::Matrix r1Swapped = la::RowSwap(r,0,1);
    la::PrintMatrix(r1Swapped);

    // -------- RowReplacement --------
    std::cout << "-------- RowReplacement --------\n\n";

    std::cout << "Let's do row replacement on row 2, using row 1 as target...\n";
    la::Matrix rr = la::RowReplacement(r,1,0);
    std::cout << "New row-equivalent matric becomes -> \n\n";
    la::PrintMatrix(rr);

    // -------- RowEchelonForm --------
    std::cout << "-------- RowEchelonForm --------\n\n";
    
    la::Vector ref1(std::vector<double>({0,2,-3,1}));
    la::Vector ref2(std::vector<double>({4,-2,1,7}));
    la::Vector ref3(std::vector<double>({-2,1,5,-8}));
    la::Matrix findref({ref1,ref2,ref3});
    std::cout << "Find REF for the below matrix...\n\n";
    la::PrintMatrix(findref);
    la::Matrix ref = la::RowEchelonForm(findref);
    std::cout << "\n\nREF is...\n\n";
    la::PrintMatrix(ref);

    // -------- ReducedRowEchelonForm --------
    std::cout << "\n-------- ReducedRowEchelonForm --------\n\n";

    std::cout << "RREF for above Ref matrix is...\n\n";
    la::Matrix rref = la::ReducedRowEchelonForm(ref);
    la::PrintMatrix(rref);

    std::cout << std::endl;
    return 0;
}