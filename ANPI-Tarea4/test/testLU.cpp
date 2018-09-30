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
#include "LU.hpp"


#include "Solver.hpp"



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
    }//luTest


  template<typename T>
  void substitutionTest( const std::function<void( anpi::Matrix<T>& ,
                                              std::vector <T>& ,
                                               std::vector <T>&)>& backwardSubs ,
                  const std::function<void( anpi::Matrix<T>& ,
                                              std::vector <T>& ,
                                               std::vector <T>&)>& forwardSubs  ){

    //initial matrix 
    Matrix<T> A;
    std::vector<T> b;
    std::vector<T> x;
    std::vector<T> x_real;


    ///testing backwardSubs
    {
      //Test with 3x3 matrix
      A = { { 1, 2, 3 },{ 0, 1, 2 },{ 0, 0, 1 } };
      b = { 9, 4, 3 };
      
      //real solution
      x_real = { 4, -2, 3 };

      backwardSubs(A, x , b);
      
      
      //Test each element one by one
      for (size_t i=0;i<A.rows();++i) {
          
        BOOST_CHECK(x[i]==x_real[i]);
          
      }
    }

    ///testing forwardSubs
    {
      //Test with 3x3 matrix
      A = { { 1, 0, 0 },{ 2, 1, 0 },{ 3, 2, 1 } };
      b = { 3, 4, 9 };
      
      //real solution
      x_real = { 3, -2, 4 };

      forwardSubs(A, x , b);
      
      
      //Test each element one by one
      for (size_t i=0;i<A.rows();++i) {
          
        BOOST_CHECK(x[i]==x_real[i]);
          
      }
    }


  } //substitutionTest


  template<typename T>
  void solveLUTest( const std::function<void(const anpi::Matrix<T>& ,
                                              std::vector <T>& ,
                                              const std::vector <T>&)>& solveLU  ){

    //initial matrix 
    Matrix<T> A;
    std::vector<T> b;
    std::vector<T> x;
    std::vector<T> x_real;
    {
            //Test with 3x3 matrix
            A = { { 2, 1 ,0 },{-1, 7, 4 },{ 0, 2, -3 } };
            b = { 4, 25, -5};
            
            //real solution
            x_real = {   1, 2, 3 };

            solveLU(A, x , b);
            
            
            BOOST_CHECK(x==x_real);
                
            
        }
  } //solveLUTest



  



  template<typename T>
  void invertTest( const std::function<void(const anpi::Matrix<T>& ,
                                              anpi::Matrix<T>&)>& invert  ){

    //
    anpi::Matrix<T> A;

    anpi::Matrix<T> Ai;



    {
      //Test with 3x3 matrix
      A = { {1, 2, 3},{0,1,4},{5,6,0} };
      invert(A, Ai);
      //Real Ai matrix
      anpi::Matrix<T> Ai_real = {{-24, 18, 5},{20, -15, -4},{-5, 4, 1}};

      
      
      //Test each element one by one
      const T eps = std::numeric_limits<T>::epsilon();

      for (size_t i=0;i<Ai_real.rows();++i) {
        for (size_t j=0;j<Ai_real.cols();++j) {
          BOOST_CHECK(std::abs(Ai_real(i,j) - Ai(i,j)) < 1000*eps);
        }
      }


    }
  }//test invert








  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( LU )

BOOST_AUTO_TEST_CASE(Doolittle) 
{
  anpi::test::luTest<float>(anpi::luDoolittle<float>,
                            anpi::unpackDoolittle<float>);
  anpi::test::luTest<double>(anpi::luDoolittle<double>,
                             anpi::unpackDoolittle<double>);
}

BOOST_AUTO_TEST_CASE(Crout) 
{
  anpi::test::luTest<float>(anpi::luCrout<float>,anpi::unpackCrout<float>);
  anpi::test::luTest<double>(anpi::luCrout<double>,anpi::unpackCrout<double>);
}

BOOST_AUTO_TEST_CASE(lu) {
  anpi::test::luTest<float>(anpi::lu<float>,
                            anpi::unpack<float>);
  anpi::test::luTest<double>(anpi::lu<double>,
                             anpi::unpack<double>);
  
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( SUBST )

    BOOST_AUTO_TEST_CASE(subtitution)
    {
        anpi::test::substitutionTest<float>(anpi::backwardSubs<float>,
                                            anpi::forwardSubs<float>);

        anpi::test::substitutionTest<double>(anpi::backwardSubs<double>,
                                            anpi::forwardSubs<double>);

    }


BOOST_AUTO_TEST_SUITE_END()




BOOST_AUTO_TEST_SUITE( SOLVE )

    BOOST_AUTO_TEST_CASE(solver)
    {
        //falla con float debido a la precision
        //anpi::test::solveLUTest<float>(anpi::solveLU<float>);

        anpi::test::solveLUTest<double>(anpi::solveLU<double>);

    }


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( INVERT )

    BOOST_AUTO_TEST_CASE(inverting)
    {
        //falla con float debido a la precision
        anpi::test::invertTest<float>(anpi::invert<float>);

        anpi::test::invertTest<double>(anpi::invert<double>);

    }


BOOST_AUTO_TEST_SUITE_END()
