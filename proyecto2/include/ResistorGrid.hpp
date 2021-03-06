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
#include "Solver.hpp"
#include "Utilities.hpp"


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

    struct Node {

        /// Row of the node
        float row_;
        /// Column of the node
        float col_;

        ///Constructors

        Node(){
        this->row_ = 0;
        this->col_ = 0;
    }

        Node(float row, float col) {
          this->row_ = row;
          this->col_ = col;
        }
    };



  class ResistorGrid
  {
  private:
    /// Matrix of the current equation system
    Matrix<float> A_;
    /// Vector of the current equation system
    std::vector<float> b_;
    /// Vector of the current equation system
    std::vector<float> c_;
    /// Raw map data
    Matrix<float> rawMap_;
    /// Min value of a resistor
    const size_t minResistor = 1;
    /// Max value of a resistor
    const size_t maxResistor = 1000000;
    /// cv::matrix of the image
    cv::Mat_<float> rawMapCV_;

  public:
    //movement vectors of the particle
    std::vector<float> posX;
    std::vector<float> posY;
    /// Matrix of the X vectors of electric field
    Matrix<float> X_;
    /// Matrix of the Y vectors of electric field
    Matrix<float> Y_;


    /// ... constructors and other methods

    /// Constructor
    
    /**
     * Constructor that initialize a matrix 'A' and a vector 'b' with the size necesary once the rawMap_ is loaded
     */
    
    ResistorGrid() {}
    
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
  
      const size_t totalVariables = 2*rows*cols-(rows+cols);
  
      A_.allocate(totalVariables,totalVariables);
      A_.fill(0.0f);
  
      b_.resize(totalVariables);
    }

    /**
     * Function that returns the index of the resistor at the edges of the grid
     * @param row1 row of the node left of the resistor
     * @param col1 column of the node left of the resistor
     * @param row2 row of the node right of the resistor
     * @param col2 column of the node right of the resistor
     * @return
     */
    std::size_t edgeResistor(std::size_t row1,
                               std::size_t col1,
                               std::size_t row2,
                               std::size_t col2) {
      return (col2-1) + (col2+1)*row2 + col2*row1;
    }

    /**
     * Convert a pair of nodes coordinates to an index, unique to each resistor
     * @param row1 row of the left node
     * @param col1 column of the left node
     * @param row2 row of the right column
     * @param col2 column of the right column
     * @return index of the resistor
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
            idmax = edgeResistor(row1, n, row2, n);
          } else {               //Horizontal case
            idmax = edgeResistor(row1, (n - 1), row2, n);
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
     * Convert an index to the pair of nodes coordinates of the resistor with that index
     * @param idx index of the resistor
     * @return struct indexPair with the pair of nodes before and after the resistor located in the rawMap_
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
     * Return the value of the resistor depending on rawMap
     * High value in case there is an obstacle
     * Low value if the space is blank
     * @param idx id of the resistor to calculate the value
     * @return value of the resistor
     */
    size_t resistorValue(const std::size_t idx) {
      
      const indexPair nodesResistor = indexToNodes(idx);
      size_t resistorOhm;
  
      if(rawMap_(nodesResistor.row1,nodesResistor.col1) == 0) {
        resistorOhm = ResistorGrid::maxResistor;
      } else if(rawMap_(nodesResistor.row2,nodesResistor.col2) == 0) {
        resistorOhm = ResistorGrid::maxResistor;
      } else {
        resistorOhm = minResistor;
      }
  
      return resistorOhm;
    }

    /**
     * Construct the grid from the given file
     * @return true if successful or false otherwise
     */
    bool build(const std::string &filename) {

      // Build the name of the image in the data path
      std::string mapPath = std::string( ANPI_DATA_PATH ) + "/" + filename;

      // Read the image using the OpenCV
      cv::imread(mapPath.c_str(),
                 CV_LOAD_IMAGE_GRAYSCALE).convertTo(rawMapCV_,CV_32FC1);
      rawMapCV_ /= 255.0f; // normalize image range to 0 .. 255

      if(rawMapCV_.cols == 0 || rawMapCV_.rows == 0 || rawMapCV_.data == NULL) {
        throw anpi::Exception("Problem creating the map");
      }

      // Convert the OpenCV matrix into an anpi matrix
      // We have to use the std::allocator to avoid an exact stride
      anpi::Matrix<float,std::allocator<float> > amapTmp(static_cast<const size_t>(rawMapCV_.rows),
                                                         static_cast<const size_t>(rawMapCV_.cols),
                                                         rawMapCV_.ptr<float>());
      // And transform it to a SIMD-enabled matrix
      anpi::Matrix<float> amap(amapTmp);

      rawMap_ = amap;

      // Initialize the other data
      const size_t totalVariables = 2*rawMap_.rows()*rawMap_.cols()-(rawMap_.rows()+rawMap_.cols());
      A_.allocate(totalVariables,totalVariables);
      A_.fill(static_cast<float>(0));

      X_.allocate(rawMap_.rows(),rawMap_.cols());
      X_.fill(static_cast<float >(0));

      Y_.allocate(rawMap_.rows(),rawMap_.cols());
      Y_.fill(static_cast<float >(0));

      //rawMap has finished building at this point

      return true;
    }

    /**
     * Compute the internal data to navigate between the given nodes using the max current to reach the end node
     * @param nodes pair with the beginning node and the end node
     * @return true if the navigation was a success
     */
    bool navigateCurrent(const indexPair& nodes) {
      generateA_(nodes);
      anpi::solveLU(A_,c_,b_);
      pathFinderMaxCurrent(nodes.row1,nodes.col1,nodes.row2,nodes.col2);
      
      return true;
    }
    
    /**
     * Compute the internal data to navigate between the given nodes using the electric field to reach the end
     * @param nodes pair with the beginning node and the end node
     * @param alpha the size of the step to take
     * @return true if the navigation was a success
     */
    bool navigateField(const indexPair& nodes, float alpha) {
      generateA_(nodes);
      anpi::solveLU(A_,c_,b_);
      
      pathFinderElectricField(nodes, alpha);
      return true;
    }

    /**
     * This method compare two Nodes and returns true if these Nodes are equals, otherwise returns false
     */
    bool compareNodes(Node node1,
                      Node node2){

      return (node1.row_ == node2.row_) && (node1.col_ == node2.col_);
    }

    /**
     * Generates the matrix of the equations with the nodes and grid equations
     * @param nodes the beginning and end nodes
     * @return bool if success
     */
    bool generateA_(const indexPair& nodes) {

      size_t numEquation = 0;
      bool equationEliminated = false;
      bool omitEquation = false;

      if((nodes.row1 > rawMap_.rows() || nodes.col1 > rawMap_.cols()) ||
         (nodes.row2 > rawMap_.rows() || nodes.col2 > rawMap_.cols())) {
        throw anpi::Exception("Invalid nodes to connect the source (nodes doesn't exist)");
      }

      if(rawMap_(nodes.row1,nodes.col1) == 0 || rawMap_(nodes.row2,nodes.col2) == 0) {
        throw anpi::Exception("Invalid nodes to connect the source (nodes connected to object)");
      }

      for(size_t i = 0; i < rawMap_.rows(); ++i) { // Careful with entering an equation when A_ is full (it shouldn't)
        for(size_t j = 0; j < rawMap_.cols(); ++j) {
          { //Nodes equations
            size_t idx;

            size_t iMinus1 = i-1;
            size_t jMinus1 = j-1;
            size_t iPlus1 = i+1;
            size_t jPlus1 = j+1;

            if(nodes.row1 == i && nodes.col1 == j) {
              b_.push_back(static_cast<float>(1));
            } else if(nodes.row2 == i && nodes.col2 == j) {
              b_.push_back(static_cast<float>(-1));
            } else {
              b_.push_back(static_cast<float>(0));
              if(!(equationEliminated)) {
                if((i==0 && j==0) ||
                   (i==0 && j==rawMap_.cols()-1) ||
                   (i==rawMap_.rows()-1 && j==0) ||
                   (i==rawMap_.rows()-1 && j==rawMap_.cols()-1)) {
                  equationEliminated = true;
                  omitEquation = true;
                  b_.pop_back();
                }
              }
            }

            if(iMinus1 < rawMap_.rows()) { // Node has upper resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(iMinus1,j,i,j);
                A_(numEquation,idx) = static_cast<float>(1);
              }
            }
            if(jMinus1 < rawMap_.cols()) { // Node has left resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(i,jMinus1,i,j);
                A_(numEquation,idx) = static_cast<float>(1);
              }
            }
            if(iPlus1 < rawMap_.rows()) { // Node has lower resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(i,j,iPlus1,j);
                A_(numEquation,idx) = static_cast<float>(-1);
              }
            }
            if(jPlus1 < rawMap_.cols()) { // Node has right resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(i,j,i,jPlus1);
                A_(numEquation,idx) = static_cast<float>(-1);
              }
            }

            if (!(omitEquation)) {
              ++numEquation;
            } else {
              omitEquation = false;
            }
          }

          {//Grid
            if(i<rawMap_.rows()-1 && j<rawMap_.cols()-1) { // Don't count the last border nodes where there's no grid
              size_t idx1 = nodesToIndex(i,j,i,j+1); // Positive current
              size_t idx2 = nodesToIndex(i,j+1,i+1,j+1); // Positive current
              size_t idx3 = nodesToIndex(i+1,j,i+1,j+1); // Negative current
              size_t idx4 = nodesToIndex(i,j,i+1,j); // Negative current

              A_(numEquation,idx1) = static_cast<float>(resistorValue(idx1));
              A_(numEquation,idx2) = static_cast<float>(resistorValue(idx2));
              A_(numEquation,idx3) = (static_cast<float>(resistorValue(idx3)))*-1;
              A_(numEquation,idx4) = (static_cast<float>(resistorValue(idx4)))*-1;

              b_.push_back(static_cast<float>(0));

              ++numEquation;
            }
          }
        }
      }

      // A_ and b_ have finished building at this point      
      
      return true;
    }

    /**
     * This method finds the path from the initial node to final node
     */
    void pathFinderMaxCurrent(size_t rowInitial,
                    size_t colInitial,
                    size_t rowFinal,
                    size_t colFinal) {

      Node previousNode(rowInitial, colInitial); //Previous Node
      Node currentNode(rowInitial, colInitial);      //Node where current's path begins and after be the current node.
      Node finalNode(rowFinal, colFinal);            //Node where current's path ends.
      cv::Mat_<float> map(rawMapCV_);
      map.at<float>((int) currentNode.row_, (int) currentNode.col_) = 0.0f; //Paint in black the current node;

      while (!(compareNodes(currentNode, finalNode))) {   //The cycle ends until the nodes are equal.
        std::vector<size_t> current_Values;        //This vector saves the index of currents of each resistor in each iteration.

        if ((currentNode.row_ == 0) &&
            (currentNode.col_ == 0)) {                                     //if current Node is (0,0)
          if (compareNodes(currentNode,
                           previousNode)) {                                              //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == -1) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      0)) {        //if current-previousNode = (-1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else {                                                                  //if current-previousNode = (0,-1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          }
        } else if ((currentNode.row_ == 0) &&
                   (currentNode.col_ == rawMap_.cols() - 1)) {                 //if current Node is (0,m)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      1)) {        //if current-previousNode = (0,1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else {                                                                   //if current-previousNode = (-1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          }
        } else if ((currentNode.row_ == rawMap_.rows() - 1) &&
                   (currentNode.col_ == 0)) {                 //if current Node is (n,0)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 1) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      0)) {        //if current-previousNode = (1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else {                                                                   //if current-previousNode = (0,-1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          }
        } else if ((currentNode.row_ == rawMap_.rows() - 1) &&
                   (currentNode.col_ == rawMap_.cols() - 1)) {  //if current Node is (n,m)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      1)) {        //if current-previousNode = (0,1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          } else {                                                                   //if current-previousNode = (1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          }
        } else if (currentNode.row_ ==
                   0) {                                                             //if current Node is (0,_)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      1)) {        //if current-previousNode = (0,1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      -1)) {   //if current-previousNode = (0,-1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else {                                                                   //if current-previousNode = (-1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          }
        } else if (currentNode.col_ ==
                   0) {                                                             //if current Node is (_,0)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 1) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      0)) {        //if current-previousNode = (1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      -1)) {  //if current-previousNode = (0,-1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else {                                                                   //if current-previousNode = (-1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          }
        } else if (currentNode.row_ ==
                   rawMap_.rows() - 1) {                                              //if current Node is (n,_)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      1)) {        //if current-previousNode = (0,1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 1) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      0)) {   //if current-previousNode = (1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
          } else {                                                                   //if current-previousNode = (0,-1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          }
        } else if (currentNode.col_ ==
                   rawMap_.cols() - 1) {                                              //if current Node is (_,m)
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      1)) {        //if current-previousNode = (0,1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 1) &&
                     ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                      0)) {   //if current-previousNode = (1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else {                                                                  //if current-previousNode = (-1,0)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
          }
        } else {                                                                    //if currentNode is in any other place in the matrix
          if (compareNodes(currentNode,
                           previousNode)) {                                             //if current Node is the initial node
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                 (size_t) currentNode.col_ + 1));
            current_Values.push_back(
                    nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                 (size_t) currentNode.col_));
          } else {
            if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                 1)) {        //if current-previousNode = (0,1)
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_ + 1));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                   (size_t) currentNode.col_));
            } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 1) &&
                       ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                        0)) {        //if current-previousNode = (1,0)
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_ + 1));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                   (size_t) currentNode.col_));
            } else if (((((int) currentNode.row_) - ((int) previousNode.row_)) == 0) &&
                       ((((int) currentNode.col_) - ((int) previousNode.col_)) ==
                        -1)) {       //if current-previousNode = (0,-1)
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_ + 1,
                                   (size_t) currentNode.col_));
            } else {                                                                        //if current-previousNode = (-1,0)
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_ - 1, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_ - 1, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_));
              current_Values.push_back(
                      nodesToIndex((size_t) currentNode.row_, (size_t) currentNode.col_, (size_t) currentNode.row_,
                                   (size_t) currentNode.col_ + 1));
            }
          }
        }
        previousNode = currentNode;                                    //Assign the previous node like currentNode
        size_t biggestVal = current_Values[0];                                             //The index used for extract the biggest value of current
        for (size_t i = 1;
             i < current_Values.size(); ++i) {           //Used for find the biggest value of current in current_Values
          if (std::abs(c_[current_Values[i - 1]]) < std::abs(c_[current_Values[i]])) {
            biggestVal = current_Values[i];
          } else {
            continue;
          }
        }
        indexPair resistor = indexToNodes(biggestVal);                 //Extract the points where is the biggest current
        //Verifying which of the both points is the after point.
        if ((currentNode.row_ == resistor.row1) && (currentNode.col_ == resistor.col1)) {
          currentNode.row_ = resistor.row2;
          currentNode.col_ = resistor.col2;
        } else {
          currentNode.row_ = resistor.row1;
          currentNode.col_ = resistor.col1;
        }
        map.at<float>((int) currentNode.row_, (int) currentNode.col_) = 0.0f; //Paint in black the current node;
      }
      cv::namedWindow("Path", CV_WINDOW_NORMAL | CV_GUI_EXPANDED);
      cv::imshow("Path", map);
      cv::waitKey();
    }

    /**
     * Bilinear interpolation formula
     * @return the new value
     */
    template <typename T>
    T bilinearInterpolationAux(T f11, T f12, T f21, T f22, T x1, T x2, T y1, T y2, T xi, T yi){

      T termino0 = (xi - x2)/(x1-x2);
      T termino1 = (yi-y2)/(y1-y2);
      T termino2 = (xi-x1)/(x2-x1);
      T termino3 = (yi-y1)/(y2-y1);

      return (termino0*termino1*f11+
              termino2*termino1*f21+
              termino0*termino3*f12+
              termino2*termino3*f22);

    }

    /**
     * Bilinear interpolation using two point between the nodes, which it rounds
     * @return the value of the x vector and y vector after bilineal
     */
    template <typename T>
    Node bilinearInterpolation(T xi, T yi){

      T x1, y1, x2, y2;
      x2 = abs(static_cast<int>(ceil(xi)));
      x1 = abs(static_cast<int>(floor(xi)));
      y2 = abs(static_cast<int>(ceil(yi)));
      y1 = abs(static_cast<int>(floor(yi)));

      ///Extract the values from the X components matrix
      T f11 = X_(x1,y1);
      T f12 = X_(x1,y2);
      T f21 = X_(x2,y1);
      T f22 = X_(x2,y2);
      ///Extract the values from the Y components matrix
      T f11s = Y_(x1,y1);
      T f12s = Y_(x1,y2);
      T f21s = Y_(x2,y1);
      T f22s = Y_(x2,y2);
      
      T row = bilinearInterpolationAux(f11, f12, f21, f22, x1, x2, y1, y2, xi, yi);
      T col = bilinearInterpolationAux(f11s, f12s, f21s, f22s, x1, x2, y1, y2, xi, yi);
      Node nodo(row, col);
      return nodo;

    }

    /**
     * Calculates the components X and Y of the current
     * @return true if success
     */
    bool calculateCurrentComponents(){

      auto minMaxElement = std::minmax_element(c_.begin(), c_.end()); // Must consider negative values of current
      float normalizeFactor = (std::abs(*minMaxElement.first) > std::abs(*minMaxElement.second)) ? // Assign the greatest value
          std::abs(*minMaxElement.first) : std::abs(*minMaxElement.second);                        // No matter the sign

      for(size_t i = 0; i < rawMap_.rows(); ++i) { // Careful with entering an equation when A_ is full (it shouldn't)
        for (size_t j = 0; j < rawMap_.cols(); ++j) {
          size_t idx;

          size_t iMinus1 = i-1;
          size_t jMinus1 = j-1;
          size_t iPlus1 = i+1;
          size_t jPlus1 = j+1;

          if(iMinus1 < rawMap_.rows()) { // Node has upper resistor
            idx = nodesToIndex(iMinus1,j,i,j);
            Y_(i,j) += c_.at(idx)*-1;
          }
          if(jMinus1 < rawMap_.cols()) { // Node has left resistor
            idx = nodesToIndex(i,jMinus1,i,j);
            X_(i,j) += c_.at(idx)*-1;
          }
          if(iPlus1 < rawMap_.rows()) { // Node has lower resistor
            idx = nodesToIndex(i,j,iPlus1,j);
            Y_(i,j) += c_.at(idx)*-1;
          }
          if(jPlus1 < rawMap_.cols()) { // Node has right resistor
            idx = nodesToIndex(i,j,i,jPlus1);
            X_(i,j) += c_.at(idx)*-1;
          }

          //Normalize the current
          //The factor is multiplied by 2 for the worse case
          //When adding vector, the value can be greater than the max current, that's why it's multiply by 2
          Y_(i,j) /= 2*normalizeFactor;
          X_(i,j) /= 2*normalizeFactor;


        }
      }

      return true;
    }

    /**
     * Finds the path from the beginning node and the end node using the electric field
     * @param nodes nodes of beginning and end of the current source
     * @param alpha the step to take between interpolations
     * @return
     */
    template <typename T>
    bool pathFinderElectricField(const indexPair& nodes, T alpha) {

      calculateCurrentComponents();

      T precision = alpha/10;

      Node pi;
      pi.col_ = nodes.col1 + alpha*X_(nodes.row1,nodes.col1);
      pi.row_ = nodes.row1 + alpha*Y_(nodes.row1,nodes.col1);
      Node velocity;
      
      posX.push_back(pi.col_);
      posY.push_back(pi.row_);

      while (abs(T(nodes.row2-pi.row_)) > precision && abs(T(nodes.col2-pi.col_)) > precision) {
        velocity = bilinearInterpolation(pi.col_,pi.row_);
        pi.col_ = pi.col_ + alpha*velocity.col_;
        pi.row_ = pi.row_ + alpha*velocity.row_;
        
        posX.push_back(pi.col_);
        posY.push_back(pi.row_);        
      }

      
      return true;
    }
  };
}


#endif //PROYECTO2_RESISTORGRID_HPP
