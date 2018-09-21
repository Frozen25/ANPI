
#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include "Utilities.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"
#include <iostream>


#ifndef ANPI_LU_HPP
#define ANPI_LU_HPP





namespace anpi{

    template <typename T>
    inline void lu(const Matrix<T>& A,
                   Matrix<T>& LU,
                   std::vector<size_t>& permut){
        luDoolittle(A, LU, permut);
    }
    template <typename T>
    inline void unpack(const Matrix<T>& LU,
                       Matrix<T>& L,
                       Matrix<T>& U){
        unpackDoolittle(LU, L, U);
    }
}







#endif //ANPI_LU_HPP
