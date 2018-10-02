//
// Created by ubuntu on 01/10/18.
//

#ifndef PROYECTO2_RESISTORGRID_HPP
#define PROYECTO2_RESISTORGRID_HPP

#include <cstdlib>
#include <iostream>

#include <AnpiConfig.hpp>

#include <string>

#include <opencv2/core.hpp>    // For cv::Mat
#include <opencv2/highgui.hpp> // For cv::imread/imshow

#include <Matrix.hpp>
#include <Exception.hpp>


namespace anpi{
  /// Pack a pair of indices of the nodes of a resistor
  struct indexPair {
    /// Row of the first node
    std::size_t row1;
    /// Column of the first node
    std::size_t col1;
    /// Row of the second node
    std::size_t row2;
    /// Column of the second node
    std::size_t col2;
  };

  class ResistorGrid
  {
  private:
    /// Matrix of the currect equation system
    Matrix<float> A_;
    /// Vector of the current equation system
    std::vector<float> b_;

    /// Raw map data
    Matrix<float> rawMap_;

  public:
    /// ... constructors and other methods

    /// Constructor
    ResistorGrid();

    /// Destructor
    ~ResistorGrid();

    /**
     * Convert a pair of nodes coordinates to an index
     */
    std::size_t nodesToIndex(std::size_t row1,
                             std::size_t col1,
                             std::size_t row2,
                             std::size_t col2);

    /**
     * Convert an index to the pair of nodes coordinates
     */
    indexPair indexToNodes(std::size_t idx); //To implement

    /**
     * Construct the grid from the given file
     * @return true if successful or false otherwise
     */
    bool build(std::string filename);

    /**
     * Compute the internal data to navigate between the given nodes
     */
    bool navigate(const indexPair& nodes);
  };
}


#endif //PROYECTO2_RESISTORGRID_HPP
