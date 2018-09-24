/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"

BOOST_AUTO_TEST_SUITE( Matrix )

/// Benchmark for subition operations
template<typename T>
class benchSub {
protected:
  /// Maximum allowed size for the square matrices
  const size_t _maxSize;

  /// A large matrix holding 
  anpi::Matrix<T> _data;

  /// State of the benchmarked evaluation
  anpi::Matrix<T> _a;
  anpi::Matrix<T> _b;
  anpi::Matrix<T> _c;
public:
  /// Construct
  benchSub(const size_t maxSize)
    : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

    size_t idx=0;
    for (size_t r=0;r<_maxSize;++r) {
      for (size_t c=0;c<_maxSize;++c) {
        _data(r,c)=idx++;
      }
    }
  }

  /// Prepare the evaluation of given size
  void prepare(const size_t size) {
    assert (size<=this->_maxSize);
    this->_a=std::move(anpi::Matrix<T>(size,size,_data.data()));
    this->_b=this->_a;
  }
};

/// Provide the evaluation method for in-place subition 
template<typename T>
class benchSubInPlaceFallback : public benchSub<T> {
public:
  /// Constructor
  benchSubInPlaceFallback(const size_t n) : benchSub<T>(n) { }
  
  // Evaluate sub in-place
  inline void eval() {
    anpi::fallback::sub(this->_a,this->_b);
  }
};

/// Provide the evaluation method for on-copy subition 
template<typename T>
class benchSubOnCopyFallback : public benchSub<T> {
public:
  /// Constructor
  benchSubOnCopyFallback(const size_t n) : benchSub<T>(n) { }
  
  // Evaluate sub on-copy
  inline void eval() {
    anpi::fallback::sub(this->_a,this->_b,this->_c);
  }
};

/// Provide the evaluation method for in-place subition 
template<typename T>
class benchSubInPlaceSIMD : public benchSub<T> {
public:
  /// Constructor
  benchSubInPlaceSIMD(const size_t n) : benchSub<T>(n) { }
  
  // Evaluate sub in-place
  inline void eval() {
    anpi::simd::sub(this->_a,this->_b);
  }
};

/// Provide the evaluation method for on-copy subition 
template<typename T>
class benchSubOnCopySIMD : public benchSub<T> {
public:
  /// Constructor
  benchSubOnCopySIMD(const size_t n) : benchSub<T>(n) { }
  
  // Evaluate sub on-copy
  inline void eval() {
    anpi::simd::sub(this->_a,this->_b,this->_c);
  }
};

/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE( sub ) {

  std::vector<size_t> sizes = {  24,  32,  48,  64,
                                 96, 128, 192, 256,
                                384, 512, 768,1024,
                               1536,2048,3072,4096};

  const size_t n=sizes.back();
  const size_t repetitions=100;
  std::vector<anpi::benchmark::measurement> times;

  {
    benchSubOnCopyFallback<float>  baoc(n);

    // Measure on-copy sub
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("sub_on_copy_float_fb.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (float) fallback","r");
  }

  {
    benchSubOnCopySIMD<float>  baoc(n);

    // Measure on-copy sub
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("sub_on_copy_float_simd.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (float) simd","g");
  }
  
  {
    benchSubInPlaceFallback<float> baip(n);

    // Measure in place sub
    ANPI_BENCHMARK(sizes,repetitions,times,baip);

    ::anpi::benchmark::write("sub_in_place_float_fb.txt",times);
    ::anpi::benchmark::plotRange(times,"In-place (float) fallback","b");
  }

  {
    benchSubInPlaceSIMD<float> baip(n);

    // Measure in place sub
    ANPI_BENCHMARK(sizes,repetitions,times,baip);

    ::anpi::benchmark::write("sub_in_place_float_simd.txt",times);
    ::anpi::benchmark::plotRange(times,"In-place (float) simd","m");
  }
  
#if 0
  
  {
    benchSubOnCopy<double>  baoc(n);

    // Measure on-copy sub
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("sub_on_copy_double.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (double)","g");
  }

  {
    benchSubInPlace<double> baip(n);

    // Measure in place sub
    ANPI_BENCHMARK(sizes,repetitions,times,baip);

    ::anpi::benchmark::write("sub_in_place_double.txt",times);
    ::anpi::benchmark::plotRange(times,"In-place (double)","m");
  }
#endif
  
  ::anpi::benchmark::show();
}
  
BOOST_AUTO_TEST_SUITE_END()
