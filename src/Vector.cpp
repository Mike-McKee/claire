#include <utility>
#include "claire/la/Vector.h"

namespace la {

Vector::Vector(std::vector<double> v) : data_(std::move(v)) {}

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