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
#include "ResistorGrid.hpp"
#include "Allocator.hpp"

BOOST_AUTO_TEST_SUITE( ResistorGrid )


/**
 * Instantiate and test the methods of the Matrix class
 */
  BOOST_AUTO_TEST_CASE( VectorField ) {
  	anpi::ResistorGrid resistorGrid;
    resistorGrid.build("mapa6x7.png");

    anpi::indexPair CurrentSource;
    CurrentSource.row1 = 0;
    CurrentSource.col1 = 0;
    CurrentSource.row2 = 4;
    CurrentSource.col2 = 5;


    resistorGrid.navigateField(CurrentSource);
    ::anpi::benchmark::quiver(resistorGrid.X_,resistorGrid.Y_,"black");

    ::anpi::benchmark::quiver1D(resistorGrid.posX,resistorGrid.posY,"black");
    


}
BOOST_AUTO_TEST_SUITE_END()





