
#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include "Utilities.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"
#include "LUDoolittle.hpp"
#include <iostream>


#ifndef ANPI_LU_HPP
#define ANPI_LU_HPP

namespace anpi{

    /**
     * Using to choose between LU decomposition with SIMD xor LU decomposition without SIMD.
     * @tparam T template.
     * @param A A matrix.
     * @param LU LU matrix.
     * @param permut permut vector.
     */
    template <typename T>
    inline void lu(const Matrix<T>& A,
                   Matrix<T>& LU,
                   std::vector<size_t>& permut){

        #ifdef ANPI_ENABLE_SIMD
          #ifdef __AVX__
            anpi::simd::luDoolittleSIMD<T,  typename avx_traits<T>::reg_type>(A, LU, permut);
          #else
            luDoolittle(A, LU, permut);
          #endif
        #else
            luDoolittle(A, LU, permut);
        #endif
    }

    /**
     * @brief used to decompose LU matrix in L matrix and U matrix.
     * @tparam T template value.
     * @param LU LU matrix.
     * @param L L matrix.
     * @param U U matrix.
     */
    template <typename T>
    inline void unpack(const Matrix<T>& LU,
                       Matrix<T>& L,
                       Matrix<T>& U){
        anpi::unpackDoolittle(LU, L, U);
    }
}







#endif //ANPI_LU_HPP
