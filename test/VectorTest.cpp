#include <iostream>
#include "claire/la/Vector.h"

int main() {

    la::Vector v({1,2,3,4});

    std::cout << "Test non-empty: Vector Size = " << v.size() << std::endl;
    std::cout << "Test non-empty: Is Vector empty = " << std::boolalpha << v.isEmpty() << std::endl;
    
    la::Vector empty({});

    std::cout << "Test empty: Vector Size = " << empty.size() << std::endl;
    std::cout << "Test empty: Is Vector empty = " << std::boolalpha << empty.isEmpty() << std::endl;

    return 0;
}