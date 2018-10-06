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


