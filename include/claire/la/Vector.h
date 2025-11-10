#pragma once
#include <vector>
#include <cstddef> // size_t

namespace la {

class Vector
{
private:
    std::vector<double> data_;
public:
    // using value_type = double; May use for future
    using size_type = std::size_t;

    //-------- Constructor --------
    explicit Vector(std::vector<double> v);
    Vector(size_type n, double value); // Creates n-dimension vector with same value for each entry
    // ~Vector();

    //-------- Accessors --------
    const std::vector<double>& data() const;
    std::vector<double>::iterator begin() {return data_.begin();}
    std::vector<double>::iterator end() {return data_.end();}
    std::vector<double>::const_iterator begin() const { return data_.begin(); }
    std::vector<double>::const_iterator end()   const { return data_.end(); }
    // std::vector<double>::reverse_iterator rbegin() {return data_.rbegin();}
    // std::vector<double>::reverse_iterator rend() {return data_.rend();}

    //-------- Object Info --------
    size_type size() const;
    bool isEmpty() const;

    double operator[](size_type i) const { return data_[i]; }
    double& operator[](size_type i) { return data_[i]; }
};

}