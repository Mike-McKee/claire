#include <utility>
#include <vector>
#include "claire/la/Vector.h"

namespace la {

Vector::Vector(std::initializer_list<double> data) : data_(std::move(data)) {}

Vector::Vector(std::vector<double> v) : data_(std::move(v)) {}

Vector::Vector(Vector::size_type n, double value) : data_(n,value) {}

const std::vector<double>& Vector::data() const {
    return data_;
}

Vector::size_type Vector::size() const {
    return data_.size();
}

bool Vector::isEmpty() const{
    if (data_.size() == 0) {
        return true;
    }
    else {return false;}
}

}