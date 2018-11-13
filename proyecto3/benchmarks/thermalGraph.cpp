/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Andres Ramirez Quiros
 * @date   30.09.2018
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
#include "ThermalPlate.hpp"
#include "Allocator.hpp"

BOOST_AUTO_TEST_SUITE( thermalGraph )


/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE( graphthermal ) {


        ::anpi::benchmark::thermal();


    }
BOOST_AUTO_TEST_SUITE_END()





