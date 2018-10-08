/*
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#ifndef ANPI_INTRINSICSM_HPP
#define ANPI_INTRINSICSM_HPP

#include <cstdint>
#include "Intrinsics.hpp"
/*
 * Include the proper intrinsics headers for the current architecture
 */

namespace anpi{

    namespace simd{


    template<typename T,class regType>
    regType mm_mul(regType, regType); 

#ifdef __AVX__
    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_mul<double>(__m256d a,__m256d b) {
      return _mm256_mul_pd(a,b);
    }
    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_mul<float>(__m256 a,__m256 b) {
      return _mm256_mul_ps(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint64_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int64_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi64(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint32_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int32_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi32(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint16_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int16_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b);
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<uint8_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b); 
    }
    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_mul<int8_t>(__m256i a,__m256i b) {
      return _mm256_mullo_epi16(a,b); 
}
#endif

}//namespace simd

}//namespace anpi
  

#endif //defined
