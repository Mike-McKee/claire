#include <iostream>
#include "claire/la/VectorOps.h"
#include "claire/la/Vector.h"

int main() {
    //------- TEST isZeroVector() -------
    // la::Vector v({1,2,3,4});
    // la::Vector w({0,0,0});

    // std::cout << "Test v = ";
    // la::PrintVector(v);
    // std::cout << " , isZeroVector = " << std::boolalpha<< la::isZeroVector(v) << std::endl;

    // std::cout << "Test w = ";
    // la::PrintVector(w);
    // std::cout << " , isZeroVector = " << std::boolalpha<< la::isZeroVector(w) << std::endl;
              
    //------- TEST isSameDimension() -------
    la::Vector u({1,4,1,7});
    la::Vector v({1,2,3,4});
    la::Vector w({0,0,0});
    
    std::cout << "Test u,w --> " << std::boolalpha << la::isSameDimension(u,w) << std::endl;
    std::cout << "Test u,v --> " << std::boolalpha << la::isSameDimension(u,v) << std::endl;
    return 0;
}