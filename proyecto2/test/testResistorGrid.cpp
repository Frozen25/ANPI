//
// Created by crisptofer-pc on 02/10/18.
//

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <functional>
#include <cmath>

#include "ResistorGrid.hpp"

BOOST_AUTO_TEST_SUITE( INDEX )

  BOOST_AUTO_TEST_CASE(nodeToIndex) {

    anpi::ResistorGrid grid(3,3);
    size_t result = grid.nodesToIndex(0,0,0,1);
    BOOST_CHECK(result ==0);
    result = grid.nodesToIndex(1,0,1,1);
    BOOST_CHECK(result ==5);
    result = grid.nodesToIndex(0,2,1,2);
    BOOST_CHECK(result == 4);
    result = grid.nodesToIndex(1,1,2,1);
    BOOST_CHECK(result == 8);
    result = grid.nodesToIndex(2,1,2,2);
    BOOST_CHECK(result == 11);
    result = grid.nodesToIndex(1,0,2,0);
    BOOST_CHECK(result == 7);
  }

  BOOST_AUTO_TEST_CASE(indexToNode) {
  
    anpi::ResistorGrid grid(3,3);
    anpi::indexPair result = grid.indexToNodes(0);
    BOOST_CHECK(result.row1 == 0);
    BOOST_CHECK(result.col1 == 0);
    BOOST_CHECK(result.row2 == 0);
    BOOST_CHECK(result.col2 == 1);
    
    result = grid.indexToNodes(5);
    BOOST_CHECK(result.row1 == 1);
    BOOST_CHECK(result.col1 == 0);
    BOOST_CHECK(result.row2 == 1);
    BOOST_CHECK(result.col2 == 1);
    
    result = grid.indexToNodes(4);
    BOOST_CHECK(result.row1 == 0);
    BOOST_CHECK(result.col1 == 2);
    BOOST_CHECK(result.row2 == 1);
    BOOST_CHECK(result.col2 == 2);
    
    result = grid.indexToNodes(8);
    BOOST_CHECK(result.row1 == 1);
    BOOST_CHECK(result.col1 == 1);
    BOOST_CHECK(result.row2 == 2);
    BOOST_CHECK(result.col2 == 1);
    
    result = grid.indexToNodes(11);
    BOOST_CHECK(result.row1 == 2);
    BOOST_CHECK(result.col1 == 1);
    BOOST_CHECK(result.row2 == 2);
    BOOST_CHECK(result.col2 == 2);
    
    result = grid.indexToNodes(7);
    BOOST_CHECK(result.row1 == 1);
    BOOST_CHECK(result.col1 == 0);
    BOOST_CHECK(result.row2 == 2);
    BOOST_CHECK(result.col2 == 0);
    
  }

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE( NAVIGATE )

  BOOST_AUTO_TEST_CASE(buildMatrix) {
    {
      try {
        anpi::ResistorGrid resistorGrid;
        resistorGrid.build("mapa3X3.png");
        BOOST_CHECK_MESSAGE(false,"Wrong name map not catched");
      } catch(anpi::Exception &exc ) {
        BOOST_CHECK_MESSAGE(true,"Successfully catched exception (Wrong map)");
      }
    }

    {
      anpi::ResistorGrid resistorGrid;
      resistorGrid.build("mapa3x3.png");

      anpi::indexPair CurrentSource;
      CurrentSource.row1 = 0;
      CurrentSource.col1 = 0;
      CurrentSource.row2 = 2;
      CurrentSource.col2 = 2;
      resistorGrid.generateA_(CurrentSource);
      BOOST_CHECK_MESSAGE(true,"Successful run 3x3");

      CurrentSource.col2 = 100;
      try {
        resistorGrid.generateA_(CurrentSource);
        BOOST_CHECK_MESSAGE(false,"Wrong source node not catched");
      } catch(anpi::Exception &exc ) {
        BOOST_CHECK_MESSAGE(true,"Successfully catched exception (Node doesn't exist)");
      }
    }

    {
      anpi::ResistorGrid resistorGrid;
      resistorGrid.build("mapa6x7.png");

      anpi::indexPair CurrentSource;
      CurrentSource.row1 = 0;
      CurrentSource.col1 = 0;
      CurrentSource.row2 = 4;
      CurrentSource.col2 = 5;
      resistorGrid.generateA_(CurrentSource);
      BOOST_CHECK_MESSAGE(true,"Successful run 6x7");
    }
  }

  BOOST_AUTO_TEST_CASE(navigateMap) {
    anpi::ResistorGrid resistorGrid;
    resistorGrid.build("mapa.png");
    
    anpi::indexPair CurrentSource;
    CurrentSource.row1 = 0;
    CurrentSource.col1 = 0;
    CurrentSource.row2 = 48;
    CurrentSource.col2 = 48;
    
    resistorGrid.navigate(CurrentSource);
  }

BOOST_AUTO_TEST_SUITE_END()

