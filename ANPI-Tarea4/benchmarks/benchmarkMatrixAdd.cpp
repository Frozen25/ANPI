/**
 * Copyright (C) 2017 
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
#include "LUDoolittle.hpp"
#include "LUCrout.hpp"


BOOST_AUTO_TEST_SUITE( MatrixDescomposition )

/// Benchmark for addition operations
    template<typename T>
    class benchDescomposition {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _maxSize;

        /// A large matrix holding
        anpi::Matrix<T> _data;

        /// State of the benchmarked evaluation
        anpi::Matrix<T> _a;
        anpi::Matrix<T> _b;
        anpi::Matrix<T> _c;
        anpi::Matrix<T> _LU;
        std::vector<size_t> _p;
    public:
        /// Construct
        benchDescomposition(const size_t maxSize)
                : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

            size_t idx=10;
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

/// Provide the evaluation method for in-place addition
    template<typename T>
    class benchLUCrout : public benchDescomposition<T> {
    public:
        /// Constructor
        benchLUCrout(const size_t n) : benchDescomposition<T>(n) { }

        // Evaluate add in-place
        inline void eval() {
            anpi::luCrout(this->_a,this->_LU,this->_p);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class benchLUDoolittle : public benchDescomposition<T> {
    public:
        /// Constructor
        benchLUDoolittle(const size_t n) : benchDescomposition<T>(n) { }

        // Evaluate add on-copy
        inline void eval() {
            anpi::luDoolittle(this->_a,this->_LU,this->_p);
        }
    };



/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE( Descomposition ) {

        std::vector<size_t> sizes = {  24,  32,  48,  64,
                                       96, 128, 192, 256};//,
                                      // 384, 512, 768,1024};
        //1536};//,2048,3072,4096};

        const size_t n=sizes.back();
        const size_t repetitions=10;
        std::vector<anpi::benchmark::measurement> times;


        {
            benchLUDoolittle<double>  baoc(n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes,repetitions,times,baoc);

            ::anpi::benchmark::write("lu_doolittle_anpi.txt",times);
            ::anpi::benchmark::plotRange(times,"lu (double)","g");
        }
        {
            benchLUCrout<double> baip(n);

            // Measure in place add
            ANPI_BENCHMARK(sizes,repetitions,times,baip);

            ::anpi::benchmark::write("lu_Crout.txt",times);
            ::anpi::benchmark::plotRange(times," lu (double) simd","m");
        }

        ::anpi::benchmark::show();
    }

BOOST_AUTO_TEST_SUITE_END()
