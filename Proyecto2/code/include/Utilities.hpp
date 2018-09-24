  
#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include <iostream>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_UTILITIES_HPP
#define ANPI_UTILITIES_HPP

namespace anpi {

///Matrices
  template<typename T>
  void swapRows(Matrix<T>& A, size_t row1Index, size_t row2Index, size_t start){
    if(row1Index != row2Index){
      for(size_t i = start; i <  A.cols(); ++i){
        T temp = A[row1Index][i];
        A[row1Index][i] = A[row2Index][i];
        A[row2Index][i] = temp;
      }
    }
  }



  template<typename T>
  void pivot(Matrix<T>& A, size_t columnIndex, size_t columnStart, size_t rowStart, std::vector<size_t>& permut){
    //Finds the maximun element in the first column to do the pivot
    T max = A[rowStart][columnIndex];
    size_t maxI = rowStart;
    for(size_t p = rowStart + 1; p < A.rows(); ++p){
      if(std::abs(A[p][columnIndex]) > std::abs(max)){
        maxI = p;
        max = A[p][columnIndex];
      }
    }
    //Swaps the row in the A matrix and in the vector
    if(maxI != rowStart){
      std::swap(permut[rowStart], permut[maxI]);
      swapRows(A, rowStart, maxI, columnStart);
    }
  }

  //this function prints a matrix
  template<typename T>
  void matrix_show(const Matrix<T>&  m, const std::string& str="") {
      std::cout << str << "\n";
      for(size_t i = 0; i < m.rows(); i++) {
          for (size_t j = 0; j < m.cols(); j++) {
              printf(" %8.3f", m(i,j));
          }
          printf("\n");
      }
      printf("\n");
  }

  template<typename T>
  anpi::Matrix<T> identityMatrix(const size_t rows, const size_t cols) {

    anpi::Matrix<T> identity(rows, cols);
    identity.fill(T(0));

    for (size_t i = 0; i < rows and i < cols; ++i) {
      identity(i,i) = T(1);
    }

    return identity;
  }

}

#endif