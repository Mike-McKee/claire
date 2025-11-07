#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/Vector.h"

// To Compile: g++ .\VectorOpsTest.cpp -I..\include\ ..\src\Vector.cpp ..\src\VectorOps.cpp -std=c++20 -o VectorOpsTest

int main() {

    //------- TEST isZeroVector() -------
    std::cout << "------- 1. TESTING isZeroVector -------\n\n";

    la::Vector v({1,2,3,4});
    la::Vector w({0,0,0});

    std::cout << "Test v = ";
    la::PrintVector(v);
    std::cout << " , isZeroVector = " << std::boolalpha<< la::isZeroVector(v) << std::endl;

    std::cout << "Test w = ";
    la::PrintVector(w);
    std::cout << " , isZeroVector = " << std::boolalpha<< la::isZeroVector(w) << std::endl;
              
    //------- TEST isSameDimension() -------
    std::cout << "\n------- 2. TESTING isSameDimension() -------\n\n";
    
    la::Vector u({1,4,1,7});
    la::Vector t({1,2,3,4});
    la::Vector s({0,0,0});
    
    std::cout << "Test u = ";
    la::PrintVector(u);
    std::cout << " and t = ";
    la::PrintVector(t);
    std::cout << " --> " << std::boolalpha << la::isSameDimension(u,t) << std::endl;

    std::cout << "Test u = ";
    la::PrintVector(u);
    std::cout << " and s = ";
    la::PrintVector(s);
    std::cout << " --> " << std::boolalpha << la::isSameDimension(u,s) << std::endl;
    

    //------- TEST VectorAddition() -------
    std::cout << "\n------- 3. TESTING VectorAddition() -------\n\n";
    
    std::cout << "Test u = ";
    la::PrintVector(u);
    std::cout << " and t = ";
    la::PrintVector(t);
    std::cout << " --> u + t = ";
    PrintVector(la::VectorAddition(u,t));
    std::cout << std::endl;
    
    //------- TEST VectorSubtraction() -------
    std::cout << "\n------- 4. TESTING VectorAddition() -------\n\n";
    
    std::cout << "Test t = ";
    la::PrintVector(t);
    std::cout << " and u = ";
    la::PrintVector(u);
    std::cout << " --> t - u = ";
    PrintVector(la::VectorSubtraction(t,u));
    std::cout << std::endl;
    
    //------- TEST ScalarMultiplication() -------
    std::cout << "\n------- 5. TESTING ScalarMultiplication() -------\n\n";
    
    std::cout << "Test u = ";
    la::PrintVector(u);
    std::cout << ", and u * 4 = ";
    PrintVector(la::ScalarMultiplication(u,4.0));
    std::cout << std::endl;
    
    //------- TEST DotProduct() -------
    std::cout << "\n------- 6. TESTING DotProduct() -------\n\n";
    
    la::Vector v1({3,-2,1});
    la::Vector v2({1,4,-2});


    std::cout << "Test v1 = ";
    la::PrintVector(v1);
    std::cout << " and v2 = ";
    la::PrintVector(v2);
    std::cout << ", then v1 dot v2 = " << la::DotProduct(v1,v2) << "\n";

    //------- TEST CrossProduct() -------
    std::cout << "\n------- 7. TESTING CrossProduct() -------\n\n";

    la::Vector v3({2,-1,3});
    la::Vector v4({1,4,-2});

    std::cout << "Test v3 = ";
    la::PrintVector(v3);
    std::cout << " and v4 = ";
    la::PrintVector(v4);
    std::cout << ", then v3 cross v4 = ";
    la::PrintVector(la::CrossProduct(v3,v4));
    std::cout << "\n\n";

    la::Vector a({3,0,-1});
    la::Vector b({-2,4,1});

    std::cout << "Test a = ";
    la::PrintVector(a);
    std::cout << " and b = ";
    la::PrintVector(b);
    std::cout << ", then a cross b = ";
    la::PrintVector(la::CrossProduct(a,b));
    std::cout << "\n";

    //------- TEST VectorNorm() -------
    std::cout << "\n------- 8. TESTING VectorNorm() -------\n\n";
    
    std::cout << "testing v3 cross v4 = ";
    la::PrintVector(la::CrossProduct(v3,v4));
    std::cout << ", and VectorNorm(v3 cross v4) = " << la::VectorNorm(la::CrossProduct(v3,v4)) << "\n";

    std::cout << "testing a cross b = ";
    la::PrintVector(la::CrossProduct(a,b));
    std::cout << ", and VectorNorm(a cross b) = " << la::VectorNorm(la::CrossProduct(a,b)) << "\n";
 
    //------- TEST Distance() -------
    std::cout << "\n------- 9. TESTING Distance() -------\n\n";

    la::Vector d1({2,-1,3});
    la::Vector d2({-1,4,1});
    
    std::cout << "Distance between d1 = ";
    la::PrintVector(d1);
    std::cout << "and d2 = ";
    la::PrintVector(d2);
    std::cout << " is -> " << la::Distance(d1,d2) << "\n";

    //------- TEST Normalize() -------
    std::cout << "\n------- 10. TESTING Normalize() -------\n\n";

    la::Vector vn({6,-3,2});
    std::cout << "The normalized vector of vn = ";
    la::PrintVector(vn);
    std::cout << " is -> ";
    la::PrintVector(la::Normalize(vn));
    std::cout << std::endl;

    //------- TEST AngleBetween() -------
    std::cout << "\n------- 11. TESTING AngleBetween() -------\n\n";

    la::Vector a1({3,2,-1});
    la::Vector a2({1,-1,2});
    std::cout << "The angle between a1 = ";
    la::PrintVector(a1);
    std::cout << " and angle a2 = ";
    la::PrintVector(a2);
    std::cout << " is -> " << la::AngleBetween(a1,a2) << "\n";

    //------- TEST Projection() -------
    std::cout << "\n------- 12. TESTING Projection() -------\n\n";

    la::Vector p1({4,2,-4});
    la::Vector p2({2,0,-1});

    std::cout << "The projection of p1 = ";
    la::PrintVector(p1);
    std::cout << " onto p2 = ";
    la::PrintVector(p2);
    std::cout << " is -> ";
    la::PrintVector(la::Projection(p1,p2));
    std::cout << std::endl;

    std::cout << std::endl;
    return 0;
}