  
#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include <iostream>
#include <string>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_UTILITIES_HPP
#define ANPI_UTILITIES_HPP

namespace anpi {



//-----------------------------------------------
///Matrices
//-----------------------------------------------
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


  #ifdef ANPI_ENABLE_SIMD
    #ifdef __AVX__
    template<typename T,typename regType>
    void swapRowsSIMD(Matrix<T>& LU,size_t r1,size_t r2){

      //std::cout << "-\n-\n-\n" << "sizeof T: "<< sizeof(T)
      //    <<"\nsize of regType: "<< sizeof(regType) <<"-\n-\n-\n"<< std::endl;
      
      unsigned long int regSize = sizeof(regType);


      //temporary element        
      regType element;
      regType* r1ptr = reinterpret_cast<regType*>(LU[r1]);
      regType* r2ptr = reinterpret_cast<regType*>(LU[r2]);


      ///total size in bytes of a row
      long unsigned int colsXsize = (long unsigned int) (LU.cols()) * (long unsigned int)(sizeof(T));

      //size_t limit = LU.cols();
      for(unsigned long int i = 0; i < colsXsize; i+=regSize ){
        //swaping the current element block at index i, between the rows.
        element = *r1ptr;
        *r1ptr++ = *r2ptr;
        *r2ptr++ = element;
      }    
      
    }
    #endif
  #endif



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

      #ifdef ANPI_ENABLE_SIMD
        #ifdef __AVX__
          swapRowsSIMD<T, typename avx_traits<T>::reg_type >(A, rowStart, maxI );
        #else
          swapRows(A, rowStart, maxI, columnStart);
        #endif
      #else
          swapRows(A, rowStart, maxI, columnStart);
      #endif
      
      std::swap(permut[rowStart], permut[maxI]);

    }
  }

  //this function prints a matrix
  template<typename T>
  static void matrix_show(const Matrix<T>&  m, const std::string& str="") {
      std::cout << str << "\n";
      for(size_t i = 0; i < m.rows(); i++) {
          for (size_t j = 0; j < m.cols(); j++) {
              printf(" %8.15f", m(i,j));
          }
          printf("\n");
      }
      printf("\n");
  }

  //this function prints a matrix
  template<typename T>
  std::string pymat_row(const Matrix<T>&  m, size_t i) {      
      std::string pyrow = "[";
      size_t j = 0;          
      for (; j < m.cols()-1; j++) {
          pyrow += std::to_string(m(i,j)) + " , ";
          
      }
      pyrow += std::to_string(m(i,j)) + " ]";
      return pyrow;
  
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

  //function prints vector
  template<typename T>
  static void vector_show(const std::vector<T> &vect){
    size_t size = vect.size();
    for(size_t i = 0; i< size ; ++i){
        std::cout << vect[i] << " " ;
      }
    std::cout << "\n";
  }

  //this function creates a new matrix Y, which is double the size of the original inside its borders
  //  and fills it with the old matrix data, for every element in A there are 4 elements in Y
  template<typename T>
  void scale_matrix(const Matrix<T>&  A, Matrix<T>&  Y) {
  size_t rows = A.cols();
  size_t cols = A.rows();

  // creates a matrix with a size equal to the double of the amount of elements outside of the border
  Y.allocate(((rows-2)*2+2),((cols-2)*2+2));

  size_t i = 0;
  size_t j = 0;
  size_t yi = 0;
  size_t yj = 0;

  for (i = 1; i < rows-1; ++i){
    // maps the old index to the corresponding index in the new matrix
    yi = (i-1)*2+1;

    for( j = 1; j < cols-1; ++j ){
      
      // maps the old index to the corresponding index in the new matrix
      yj = (j-1)*2+1;

      //stores old values in new matrix
      Y[yi ][yj ] = A[i][j];
      Y[yi ][yj + 1] = A[i][j];
      Y[yi + 1 ][yj ] = A[i][j];
      Y[yi + 1][yj + 1] = A[i][j];
      }
    }
  }






}

#endif