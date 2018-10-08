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

    struct Node {

        /// Row of the node
        std::size_t row_;
        /// Column of the node
        std::size_t col_;


        ///Constructors

        Node(){
        this->row_ = 0;
        this->col_ = 0;
    }

        Node(size_t row, size_t col) {
          this->row_ = row;
          this->col_ = col;
        }
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
    /// Min value of a resistor
    const size_t minResistor = 1;
    /// Max value of a resistor
    const size_t maxResistor = 1000000;

  public:
    /// ... constructors and other methods

    /// Constructor
    
    /**
     * Constructor that initialize a matrix 'A' and a vector 'b' with the size necesary once the rawMap_ is loaded
     */
    ResistorGrid() {
      const size_t totalVariables = 2*rawMap_.rows()*rawMap_.cols()-(rawMap_.rows()+rawMap_.cols());
      A_.allocate(totalVariables,totalVariables);
      A_.fill(0.0f);
  
      b_.resize(totalVariables);
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
  
      try {
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
        
      } catch (const std::exception &exc){
        // catch anything thrown within try block that derives from std::exception
        std::cerr << exc.what();
        return false;
      }
    }

    /**
     * Compute the internal data to navigate between the given nodes
     */
    bool navigate(const indexPair& nodes) {
      
      size_t numEquation = 0;
      const size_t totalVariables = 2*rawMap_.rows()*rawMap_.cols()-(rawMap_.rows()+rawMap_.cols());
      bool equationEliminated = false;
      bool omitEquation = false;
      
      for(size_t i = 0; i < rawMap_.rows(); ++i) {
        for(size_t j = 0; i < rawMap_.cols(); ++j) {
          { //Nodes equations
            size_t idx;
            
            size_t iMinus1 = i-1;
            size_t jMinus1 = j-1;
            size_t iPlus1 = i+1;
            size_t jPlus1 = j+1;
  
            if(nodes.row1 == i && nodes.col1 == j) {
              b_.push_back(1);
            } else if(nodes.row2 == i && nodes.col2 == j) {
              b_.push_back(-1);
            } else {
              b_.push_back(0);
              if(!(equationEliminated)) {
                if((i==0 && j==0) ||
                   (i==0 && j==rawMap_.cols()-1) ||
                   (i==rawMap_.rows()-1 && j==0) ||
                   (i==rawMap_.rows()-1 && j==rawMap_.cols()-1)) {
                  equationEliminated = true;
                  omitEquation = true;
                }
              }
            }
  
            if(iMinus1 < totalVariables) { // Node has upper resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(iMinus1,j,i,j);
                A_(numEquation,idx) = 1;
              }
            }
            if(jMinus1 < totalVariables) { // Node has left resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(i,jMinus1,i,j);
                A_(numEquation,idx) = 1;
              }
            }
            if(iPlus1 < totalVariables) { // Node has lower resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(i,j,iPlus1,j);
                A_(numEquation,idx) = -1;
              }
            }
            if(jPlus1 < totalVariables) { // Node has right resistor
              if(!(omitEquation)) {
                idx = nodesToIndex(i,j,i,jPlus1);
                A_(numEquation,idx) = -1;
              }
            }
  
            if (!(omitEquation)) {
              ++numEquation;
            } else {
              omitEquation = false;
            }
          }
          
          {//Grid
            size_t idx1;
            size_t idx2;
            size_t idx3;
            size_t idx4;
            
          }
        }
      }
    }
    /**
     * This method compare two Nodes and returns true if these Nodes are equals, otherwise returns false
     */
    bool compareNodes(Node node1,
                      Node node2){

      return (node1.row_ == node2.row_) && (node1.col_ == node2.col_);
    }


    /**
     * This method finds the path from the initial node to final node
     */
    template <typename T>
    void pathFinder(size_t rowInitial,
                    size_t colInitial,
                    size_t rowFinal,
                    size_t colFinal,
                    const std::vector<T>& SolVector){

      Node previousNode(rowInitial, colInitial); //Previous Node
      Node currentNode(rowInitial, colInitial);      //Node where current's path begins and after be the current node.
      Node finalNode(rowFinal, colFinal);            //Node where current's path ends.
      rawMap_[currentNode.row_][currentNode.col_] = 0;  //Paint in black the current node;

      while(!(compareNodes(currentNode, finalNode))){   //The cycle ends until the nodes are equal.

        std::vector<size_t > current_Values;        //This vector saves the index of currents of each resistor in each iteration.


        if((currentNode.row_ == 0) && (currentNode.col_ == 0)){                                     //if current Node is (0,0)

          if(compareNodes(currentNode, previousNode)){                                              //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));
          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == -1) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 0)) {        //if current-previousNode = (-1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }
          else {                                                                  //if current-previousNode = (0,-1)
            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));
          }

        }

        else if((currentNode.row_ == 0) && (currentNode.col_ == rawMap_.cols()-1)){                 //if current Node is (0,m)

          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1,  currentNode.row_, currentNode.col_));

          }
          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 1)) {        //if current-previousNode = (0,1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

          }
          else{                                                                   //if current-previousNode = (-1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1,  currentNode.row_, currentNode.col_));

          }

        }

        else if((currentNode.row_ == rawMap_.rows()-1) && (currentNode.col_ == 0)){                 //if current Node is (n,0)


          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 1) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 0)) {        //if current-previousNode = (1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }
          else{                                                                   //if current-previousNode = (0,-1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));

          }

        }

        else if((currentNode.row_ == rawMap_.rows()-1) && (currentNode.col_ == rawMap_.cols()-1)){  //if current Node is (n,m)


          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 1)) {        //if current-previousNode = (0,1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));

          }
          else{                                                                   //if current-previousNode = (1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));

          }

        }

        else if(currentNode.row_ == 0){                                                             //if current Node is (0,_)

          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 1)) {        //if current-previousNode = (0,1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }
          else if (((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
                   ((((int)currentNode.col_)-((int)previousNode.col_)) == -1)){   //if current-previousNode = (0,-1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

          }
          else{                                                                   //if current-previousNode = (-1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }

        }

        else if(currentNode.col_ == 0){                                                             //if current Node is (_,0)

          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));


          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 1) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 0)) {        //if current-previousNode = (1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
                   ((((int)currentNode.col_)-((int)previousNode.col_)) == -1)) {  //if current-previousNode = (0,-1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

          }
          else{                                                                   //if current-previousNode = (-1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }

        }

        else if(currentNode.row_ == rawMap_.rows()-1){                                              //if current Node is (n,_)

          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));


          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 1)) {        //if current-previousNode = (0,1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));


          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 1) &&
                   ((((int)currentNode.col_)-((int)previousNode.col_)) == 0)) {   //if current-previousNode = (1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

          }
          else{                                                                   //if current-previousNode = (0,-1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));

          }

        }

        else if(currentNode.col_ == rawMap_.cols()-1){                                              //if current Node is (_,m)

          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));


          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
              ((((int)currentNode.col_)-((int)previousNode.col_)) == 1)) {        //if current-previousNode = (0,1)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

          }

          else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 1) &&
                   ((((int)currentNode.col_)-((int)previousNode.col_)) == 0)) {   //if current-previousNode = (1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

          }

          else {                                                                  //if current-previousNode = (-1,0)

            //Add the corresponding index of current in the vector.
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));

          }

        }
        else{                                                                    //if currentNode is in any other place in the matrix

          if(compareNodes(currentNode, previousNode)) {                                             //if current Node is the initial node

            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));
            current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));


          }
          else{

            if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
                ((((int)currentNode.col_)-((int)previousNode.col_)) == 1)) {        //if current-previousNode = (0,1)

              current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));
              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

            }

            else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 1) &&
                     ((((int)currentNode.col_)-((int)previousNode.col_)) == 0)) {        //if current-previousNode = (1,0)

              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));
              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

            }

            else if( ((((int)currentNode.row_)-((int)previousNode.row_)) == 0) &&
                     ((((int)currentNode.col_)-((int)previousNode.col_)) == -1)) {       //if current-previousNode = (0,-1)

              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
              current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_+1, currentNode.col_));

            }

            else{                                                                        //if current-previousNode = (-1,0)

              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_-1, currentNode.row_, currentNode.col_));
              current_Values.push_back(nodesToIndex(currentNode.row_-1, currentNode.col_, currentNode.row_, currentNode.col_));
              current_Values.push_back(nodesToIndex(currentNode.row_, currentNode.col_, currentNode.row_, currentNode.col_+1));

            }
          }
        }

        previousNode = currentNode;                                    //Assign the previous node like currentNode
        size_t biggestVal;                                             //The index used for extract the biggest value of current
        for (size_t i = 0; i < current_Values.size(); ++i) {           //Used for find the biggest value of current in current_Values

          if(i == 0){
            biggestVal = current_Values[i];
          }
          else{

            if(SolVector[current_Values[i-1]] < SolVector[current_Values[i]]){

              biggestVal = current_Values[i];

            }
            else{
              continue;
            }
          }
        }
        indexPair resistor = indexToNodes(biggestVal);                 //Extract the points where is the biggest current

        //Verifying which of the both points is the after point.
        if((currentNode.row_ == resistor.row1) && (currentNode.col_ == resistor.col1)){
          currentNode.row_ = resistor.row2;
          currentNode.col_ = resistor.col2;
        } else{
          currentNode.row_ = resistor.row1;
          currentNode.col_ = resistor.col1;
        }

        rawMap_[currentNode.row_][currentNode.col_] = 0;  //Paint in black the new current node;

      }
    }



  };
}


#endif //PROYECTO2_RESISTORGRID_HPP