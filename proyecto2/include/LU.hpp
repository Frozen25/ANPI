
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

        #ifdef ANPI_ENABLE_SIMD
          #ifdef __AVX__
            simd::luDoolittleSIMD<T,  typename avx_traits<T>::reg_type>(A, LU, permut);
          #else
            luDoolittle(A, LU, permut);
          #endif
        #else
            luDoolittle(A, LU, permut);
        #endif


        //luDoolittle(A, LU, permut);
    }
    template <typename T>
    inline void unpack(const Matrix<T>& LU,
                       Matrix<T>& L,
                       Matrix<T>& U){
        unpackDoolittle(LU, L, U);
    }
}







#endif //ANPI_LU_HPP
