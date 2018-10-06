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
    ResistorGrid() {
      size_t ASize = 2*rawMap_.rows()*rawMap_.cols()-(rawMap_.rows()+rawMap_.cols());
      A_.allocate(ASize,ASize);
      A_.fill(0.0f);
  
      b_.resize(ASize);
    }
    
    ResistorGrid(size_t rows, size_t cols) {
      rawMap_.allocate(rows,cols);
      rawMap_.fill(0.0f);
  
      size_t ASize = 2*rows*cols-(rows+cols);
  
      A_.allocate(ASize,ASize);
      A_.fill(0.0f);
  
      b_.resize(ASize);
    }

    std::size_t extremos(std::size_t row1,
                          std::size_t col1,
                          std::size_t row2,
                          std::size_t col2) {
      return (col2-1) + (col2+1)*row2 + col2*row1;
    }

    /**
     * Convert a pair of nodes coordinates to an index
     */

    std::size_t nodesToIndex(std::size_t row1,
                             std::size_t col1,
                             std::size_t row2,
                             std::size_t col2) {
  
      if ((std::abs(int(row2-row1)) == 1 || row2-row1 == 0) &&
          (std::abs(int(col2-col1)) == 1 || col2-col1 == 0)) {
        if (row2 == row1+1 || col2 == col1+1) {
          size_t idmax, idx;
          size_t n = rawMap_.cols()-1;
      
          if(col1 == col2) {    //Vertical case
            idmax = extremos(row1, n, row2, n);
          } else {               //Horizontal case
            idmax = extremos(row1, (n-1), row2, n);
          }
      
          idx = idmax - (n-col2);
          return idx;
        } else {
          throw anpi::Exception("Second pair must be greater than first pair");
        }
      } else {
        throw anpi::Exception("Nodes are not next to each other");
      }
    }

    /**
     * Convert an index to the pair of nodes coordinates
     */
    indexPair indexToNodes(std::size_t idx) {
  
      if (((2*rawMap_.cols()*rawMap_.rows())-
           (rawMap_.cols()+rawMap_.rows())) > idx) {
    
        indexPair result;
        size_t n = rawMap_.cols()-1;
        size_t horizontal, vertical;
    
        for (size_t i = 0;;++i) {
          vertical = nodesToIndex(i,n,i+1,n);
          horizontal = nodesToIndex(i,(n-1),i,n);
      
          if (horizontal >= idx){
            result.row1 = i;
            result.col1 = (n-1)-(horizontal-idx);
            result.row2 = i;
            result.col2 = n-(horizontal-idx);
            break;
          } else if (vertical >= idx) {
            result.row1 = i;
            result.col1 = n-(vertical-idx);
            result.row2 = i+1;
            result.col2 = n-(vertical-idx);
            break;
          }
        }
        return result;
    
      } else {
        throw anpi::Exception("Idx not in range");
      }
    }

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
