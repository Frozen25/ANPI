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

#include "ThermalPlate.hpp"
#include "Utilities.hpp"


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
  namespace test {
  

  template<typename T>
  void scalingTest( const std::function<void( anpi::Matrix<T>& ,
                                              anpi::Matrix<T>& )>& scale_matrix  ){

    //initial matrix 
    anpi::Matrix<T> A;
    anpi::Matrix<T> B;
    {
      //Test with 3x3 matrix
      A = { { 2, 3 ,0 },{-1, 2, 4 },{ 0, 2, -3 } };

      anpi::ThermalPlate thermalPlate;
      thermalPlate.scale_matrix(A, B);
      
      
      BOOST_CHECK(1==1);    
      BOOST_CHECK_MESSAGE(true,"Successful run 6x7");   
    }


    


  } //scalingTest

  
  
  } // test
}  // anpi





BOOST_AUTO_TEST_SUITE( SOLVEThermal )

    BOOST_AUTO_TEST_CASE(solverthermal)
    {
      
        //initial matrix 
      anpi::Matrix<double> A;
      anpi::Matrix<double> B;
      anpi::Matrix<double> C;
      const double eps = std::numeric_limits<double>::epsilon();
      anpi::ThermalPlate thermalPlate;


      //Test with 3x3 matrix
      A = { { 2.0f, 3.0f ,0.0f },{-1.0f, 2.0f, 4.0f },{ 0.0f, 2.0f, -3.0f } };
      
      thermalPlate.scale_matrix(A, B);
      
      anpi::Matrix<double> B_real = {
          { 0.00f, 3.00f, 3.00f, 3.00f},
          {-1.00f, 1.50f, 2.75f, 4.00f},
          {-1.00f, 1.25f, 2.50f, 4.00f},
          {-1.00f, 2.00f, 2.00f, 4.00f}};

      for(size_t i = 1; i< B.rows()-1 ;++i){
        for(size_t j = 1; j< B.cols()-1; ++j){
          BOOST_CHECK(std::abs(B_real(i,j) - B(i,j)) < 0.000001f);
        }
      }



      //Test with 5x5 matrix
      thermalPlate.scale_matrix(B, C);
      
      anpi::Matrix<double> C_real = {
           { 0.0000f, 3.0000f, 3.0000f, 3.0000f, 3.0000f, 3.0000f},
           {-1.0000f, 1.2500f, 2.1875f, 2.5000f, 3.1250f, 4.0000f},
           {-1.0000f, 0.8125f, 1.7500f, 2.3750f, 3.0000f, 4.0000f},
           {-1.0000f, 0.7500f, 1.6250f, 2.2500f, 2.9375f, 4.0000f},
           {-1.0000f, 0.8750f, 1.7500f, 2.0625f, 2.7500f, 4.0000f},
           {-1.0000f, 2.0000f, 2.0000f, 2.0000f, 2.0000f, 4.0000f}};

      for(size_t i = 1; i< C.rows()-1 ;++i){
        for(size_t j = 1; j< C.cols()-1; ++j){
          BOOST_CHECK(std::abs(C_real(i,j) - C(i,j)) < 1000*eps);
        }
      }





      //BOOST_CHECK_MESSAGE(true,"Successful run 6x7");   
    
      //anpi::test::scalingTest<double>(anpi::ThermalPlate.scale_matrix<double>);

    }


BOOST_AUTO_TEST_SUITE_END()

