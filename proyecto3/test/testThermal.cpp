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
#include <vector>




BOOST_AUTO_TEST_SUITE( ScaleThermal )

    BOOST_AUTO_TEST_CASE( thermalscaling )
    {
      
        //initial matrix 
      anpi::Matrix<float> A;
      anpi::Matrix<float> B;
      anpi::Matrix<float> C;
      const float eps = std::numeric_limits<float>::epsilon();


      std::vector<float> top = { 3.0f};
      std::vector<float> bottom = { 2.0f};
      std::vector<float> left = { -1.0f};
      std::vector<float> right = { 8.0f};
      //const float eps = std::numeric_limits<float>::epsilon();
      anpi::ThermalPlate thermalPlate(0,0,0,0,
                                    top, bottom, left, right, 0, 0,
                                    0, 0, 0);


      //Test with 3x3 matrix
      A = { { 2.0f, 3.0f ,0.0f },{-1.0f, 2.0f, 8.0f },{ 0.0f, 2.0f, -3.0f } };
      
      thermalPlate.scale_matrix(A, B);

      
      
      anpi::Matrix<float> B_real =  {
           { 0.00, 3.00, 3.00, 0.00},
           {-1.00, 1.50, 3.75, 8.00},
           {-1.00, 1.25, 3.50, 8.00},
           { 0.00, 2.00, 2.00, 0.00}};

      for(size_t i = 1; i< B.rows()-1 ;++i){
        for(size_t j = 1; j< B.cols()-1; ++j){
          BOOST_CHECK(std::abs(B_real(i,j) - B(i,j)) < 0.000001f);
        }
      }



      //Test with 5x5 matrix
      thermalPlate.scale_matrix(B, C);
      
      anpi::Matrix<float> C_real = {
           { 0.0000, 3.0000, 3.0000, 3.0000, 3.0000, 0.0000},
           {-1.0000, 1.2500, 2.4375, 3.0000, 4.6250, 8.0000},
           {-1.0000, 0.8125, 2.0000, 3.1250, 4.7500, 8.0000},
           {-1.0000, 0.7500, 1.8750, 3.0000, 4.6875, 8.0000},
           {-1.0000, 0.8750, 2.0000, 2.5625, 4.2500, 8.0000},
           { 0.0000, 2.0000, 2.0000, 2.0000, 2.0000, 0.0000}};


      for(size_t i = 1; i< C.rows()-1 ;++i){
        for(size_t j = 1; j< C.cols()-1; ++j){
          BOOST_CHECK(std::abs(C_real(i,j) - C(i,j)) < 1000*eps);
        }
      }





      //BOOST_CHECK_MESSAGE(true,"Successful run 6x7");   
    
      //anpi::test::scalingTest<float>(anpi::ThermalPlate.scale_matrix<float>);

    }


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( IterateThermal )

    BOOST_AUTO_TEST_CASE(thermaliteration)
    {
      
        //initial matrix 
      anpi::Matrix<float> A;
      anpi::Matrix<float> B;
      
      //const float eps = std::numeric_limits<float>::epsilon();
      std::vector<float> top = { 3.0f};
      std::vector<float> bottom = { 2.0f};
      std::vector<float> left = { -1.0f};
      std::vector<float> right = { 8.0f};
      //const float eps = std::numeric_limits<float>::epsilon();
      anpi::ThermalPlate thermalPlate(0,0,0,0,
                                    top, bottom, left, right, 0, 0,
                                    0, 0, 0);



      //Test with 3x3 matrix
      A = { { 2.0f, 3.0f ,0.0f },{-1.0f, 2.0f, 8.0f },{ 0.0f, 2.0f, -3.0f } };
      

      thermalPlate.calculatePlate(A, B, 1.001f);
      



         
    
      //anpi::test::scalingTest<float>(anpi::ThermalPlate.scale_matrix<float>);

    }


BOOST_AUTO_TEST_SUITE_END()





BOOST_AUTO_TEST_SUITE( SolveThermal )

    BOOST_AUTO_TEST_CASE(thermalsolution)
    {
      
        //initial matrix 
      anpi::Matrix<float> A;
      anpi::Matrix<float> B;

      /*
      
      std::vector<float> top = { 3.0f, 999.0f};
      std::vector<float> bottom = { 2.0f, 30.0f,10.0f};
      std::vector<float> left = { 10.0f, 30.0f, -10.0f, -100.0f , 220.0f};
      std::vector<float> right = { 80.0f, 50.0f, 80.0f};
      //const float eps = std::numeric_limits<float>::epsilon();
      
      */
      std::vector<float> top = { 3.0f, 220.0f};
      std::vector<float> bottom = { 0.0f , 100.0f };
      std::vector<float> left = { 10.0f, 30.0f, -10.0f, -100.0f , 220.0f};
      std::vector<float> right = { 80.0f, 50.0f, 80.0f};
      //const float eps = std::numeric_limits<float>::epsilon();
      
      
      anpi::ThermalPlate thermalPlate(0,0,0,0,
                                    top, bottom, left, right, 0, 0,
                                    0, 0, 0);
    

      //Test with 3x3 matrix
      
      //// parameters are: epsilon  iterations   saveFlag
      thermalPlate.solvePlate(2.1f, 15, 1);
      
      



         
    
      //anpi::test::scalingTest<float>(anpi::ThermalPlate.scale_matrix<float>);

    }


BOOST_AUTO_TEST_SUITE_END()