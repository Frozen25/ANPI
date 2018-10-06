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
    
    /**
     * Constructor that initialize a matrix 'A' and a vector 'b' with the size necesary once the rawMap_ is loaded
     */
    ResistorGrid() {
      size_t ASize = 2*rawMap_.rows()*rawMap_.cols()-(rawMap_.rows()+rawMap_.cols());
      A_.allocate(ASize,ASize);
      A_.fill(0.0f);
  
      b_.resize(ASize);
    }
    
    /**
     * Constructor of the class that initialize an empty map of row*cols,
     * a solution matrix 'A' of size depending of the variables
     * and a 'b' vector of size equal to 'A'
     * @param rows number of rows of rawMap_
     * @param cols number of columns of rawMap_
     */
    ResistorGrid(size_t rows, size_t cols) {
      rawMap_.allocate(rows,cols);
      rawMap_.fill(0.0f);
  
      size_t ASize = 2*rows*cols-(rows+cols);
  
      A_.allocate(ASize,ASize);
      A_.fill(0.0f);
  
      b_.resize(ASize);
    }

    /**
     * Function that returns the index of the resistance at the edges of the grid
     * @param row1 row of the node left of the resistance
     * @param col1 column of the node left of the resistance
     * @param row2 row of the node right of the resistance
     * @param col2 column of the node right of the resistance
     * @return
     */
    std::size_t edgeResistance(std::size_t row1,
                               std::size_t col1,
                               std::size_t row2,
                               std::size_t col2) {
      return (col2-1) + (col2+1)*row2 + col2*row1;
    }

    /**
     * Convert a pair of nodes coordinates to an index, unique to each resistance
     * @param row1 row of the left node
     * @param col1 column of the left node
     * @param row2 row of the right column
     * @param col2 column of the right column
     * @return index of the resistance
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
            idmax = edgeResistance(row1, n, row2, n);
          } else {               //Horizontal case
            idmax = edgeResistance(row1, (n - 1), row2, n);
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
     * Convert an index to the pair of nodes coordinates of the resistance with that index
     * @param idx index of the resistance
     * @return struct indePair with the pair of nodes before and after the resistance located in the rawMap_
     */
    indexPair indexToNodes(std::size_t idx) {
  
      if (((2*rawMap_.cols()*rawMap_.rows())-
           (rawMap_.cols()+rawMap_.rows())) > idx) {
    
        indexPair result{};
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
    bool build(const std::string &filename) {
      // Build the name of the image in the data path
      std::string mapPath = std::string( ANPI_DATA_PATH ) + filename;
  
      // Read the image using the OpenCV
      cv::Mat_<float> map;
  
      cv::imread(mapPath.c_str(),
                 CV_LOAD_IMAGE_GRAYSCALE).convertTo(map,CV_32FC1);
      map /= 255.0f; // normalize image range to 0 .. 255
  
      // Convert the OpenCV matrix into an anpi matrix
      // We have to use the std::allocator to avoid an exact stride
      anpi::Matrix<float,std::allocator<float> > amapTmp(static_cast<const size_t>(map.rows),
                                                         static_cast<const size_t>(map.cols),
                                                         map.ptr<float>());
      // And transform it to a SIMD-enabled matrix
      anpi::Matrix<float> amap(amapTmp);
  
      rawMap_ = amap;
  
      return true;
    }

    /**
     * Compute the internal data to navigate between the given nodes
     */
    bool navigate(const indexPair& nodes);
  };
}


#endif //PROYECTO2_RESISTORGRID_HPP
