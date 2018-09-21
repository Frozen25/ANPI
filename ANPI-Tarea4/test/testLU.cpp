/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>

#include "LUCrout.hpp"
#include "LUDoolittle.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
  namespace test {
    
    /// Test the given closed root finder
    template<typename T>
    void luTest(const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         std::vector<size_t>&)>& decomp,
                const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         Matrix<T>&)>& unpack) {

      // The result
      Matrix<T> LU;

      // Test if a non-square matrix is successfully detected
      {
        Matrix<T> A = {{1,7,6,4},{2,17,27,17}};
        std::vector<size_t> p;
        try {
          decomp(A,LU,p);
          BOOST_CHECK_MESSAGE(false,"Rectangular matrix not properly catched");
        }
        catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Rectangular matrix properly detected");
        }
      }

      // Test pivoting
      {
        anpi::Matrix<T> A = { {-1,-2,1,2},{ 2, 0,1,2},{-1,-1,0,1},{ 1, 1,1,1} };
        std::vector<size_t> p;
        decomp(A,LU,p);

        std::vector<size_t> gp= {1,0,3,2};
        BOOST_CHECK(gp==p);
      }
      
      // Test decomposition
      {
        // same matrix as before, but already permuted to force a
        // clean decomposition
        anpi::Matrix<T> A = { { 2, 0,1,2},{-1,-2,1,2},{ 1, 1,1,1},{-1,-1,0,1} };
        std::vector<size_t> p;
        decomp(A,LU,p);
        Matrix<T> L,U;
        unpack(LU,L,U);
        Matrix<T> Ar=L*U;
        const T eps = std::numeric_limits<T>::epsilon();
        BOOST_CHECK(Ar.rows()==A.rows());
        BOOST_CHECK(Ar.cols()==A.cols());
        for (size_t i=0;i<Ar.rows();++i) {
          for (size_t j=0;j<Ar.cols();++j) {
            BOOST_CHECK(std::abs(Ar(i,j)-A(i,j)) < eps);
          }
        }
      }

      




    }

      template<typename T>
      void unpackDoolittleTest(const std::function<void(const Matrix<T>&,
              Matrix<T>&, Matrix<T>&)>& unpack) {    //Test the unpack doolittle method

          //factorized matrix in LU form with Doolittle
        Matrix<T> LU;

        Matrix<T> L;
        Matrix<T> U;

        {
          //Test with 4x4 matrix
          LU = { {1, 2, 3, 4},{5, 6, 7, 8},{9, 10, 11, 12},{13, 14, 15, 6} };
          unpack(LU, L, U);
          //Real L matrix
          Matrix<T> Lr = {{1,0,0,0},{5,1,0,0},{9,10,1,0},{13,14,15,1}};
          //Real U matrix
          Matrix<T> Ur = {{1,2,3,4},{0,6,7,8},{0,0,11,12},{0,0,0,16}};
          //Test each element one by one
          for (size_t i=0;i<LU.rows();++i) {
            for (size_t j=0;j<LU.cols();++j) {
              BOOST_CHECK(L(i,j)==Lr(i,j));
            }
          }
        }


        {
          //Test with 3x3 matrix
          LU = { {1, 2, 3},{4, 5, 6},{7, 8, 9} };
          unpack(LU, L, U);
          //Real L matrix
          Matrix<T> Lr = {{1, 0, 0},{4, 1, 0},{7, 8, 1}};
          //Real U matrix
          Matrix<T> Ur = {{1, 2, 3},{0, 5, 6},{0, 0, 9}};
          //Test each element one by one
          for (size_t i=0;i<LU.rows();++i) {
            for (size_t j=0;j<LU.cols();++j) {
              BOOST_CHECK(L(i,j)==Lr(i,j));
            }
          }
        }
      }


  } // test
}  // anpi






BOOST_AUTO_TEST_SUITE( LU )

BOOST_AUTO_TEST_CASE(Doolittle) 
{
  anpi::test::luTest<float>(anpi::luDoolittle<float>,
                            anpi::unpackDoolittle<float>);
  anpi::test::luTest<double>(anpi::luDoolittle<double>,
                             anpi::unpackDoolittle<double>);
  anpi::test::unpackDoolittleTest<float>(anpi::unpackDoolittle<float>);

  anpi::test::unpackDoolittleTest<double>(anpi::unpackDoolittle<double>);

}

BOOST_AUTO_TEST_CASE(Crout) 
{
  anpi::test::luTest<float>(anpi::luCrout<float>,anpi::unpackCrout<float>);
  anpi::test::luTest<double>(anpi::luCrout<double>,anpi::unpackCrout<double>);
}


BOOST_AUTO_TEST_SUITE_END()
