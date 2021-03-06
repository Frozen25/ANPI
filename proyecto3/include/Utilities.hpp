  
#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_UTILITIES_HPP
#define ANPI_UTILITIES_HPP


namespace anpi {



//-----------------------------------------------
///Matrices
//-----------------------------------------------

/**
  * @brief Used to swap the rows of any matrix
  *
  * @tparam T template value
  * @param A The matrix that we want to swap yours rows
  * @param row1Index Index of the row 1
  * @param row2Index Index of the row 2
  * @param start value where we want to start.
  */
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

    /**
      * @brief Used to swap the rows of any matrix using simd instructions.
      *
      * @tparam T template value
      * @tparam regType template value of the register.
      * @param LU Matrix to we want to swap his rows.
      * @param r1 size of row 1.
      * @param r2 size of row 2.
      */
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


 /**
  * @brief used to swap some rows to reduce numerical errors.
  * @tparam T template value
  * @param A matrix to we want apply pivoting process.
  * @param columnIndex index of the column to apply the pivot.
  * @param columnStart column where to start the pivoting process.
  * @param rowStart row where to start the pivoting process.
  * @param permut permutation vector.
  */
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

  /**
   * @brief used to print the matrix
   * @tparam T template value
   * @param m matrix to print
   * @param str string to print with matrix
   */
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

  //this function prints a matrix casting the values to integer of default size 3 (max number 999)
  template<typename T>
  static void matrix_show_int(const Matrix<T>&  m, const std::string& str="", int maxsize = 3) {
      std::cout << str << "\n";
      std::stringstream ss;
      for(size_t i = 0; i < m.rows(); i++) {
          for (size_t j = 0; j < m.cols(); j++) {              
              ss.str("");
              ss << std::setw(maxsize) << std::setfill (' ') << (int)(m(i,j));
              std::cout << ss.str() << ' ';
              //printf(" %d",(int) m(i,j));
          }
          printf("\n");
      }
      printf("\n");
  }

  //this function saves a matrix casting the values to integer of default size 3 (max number 999)
  // saves the matrix in a file called matrix.txt
  template<typename T>
  static void matrix_show_file(const Matrix<T>&  m,  bool novisuals) {
      
      //std::stringstream ss;

      std::ofstream myfile;
      myfile.open ("matrix.txt");

      if (!novisuals){
        for(size_t i = 0; i < m.rows(); i++) {
            for (size_t j = 0; j < m.cols(); j++) {              
                //ss.str("");
                //ss << std::setw(maxsize) << std::setfill (' ') << (int)(m(i,j));
                //myfile << ss.str() << ' ';              
                myfile << m(i,j) << ' ';
            }
            myfile << "\n";
        }
        std::cout << "Saved matix to file: matrix.txt\n";
      }

      myfile << "\n";

      myfile.close();
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


  /**
   * @brief generate a identity matrix.
   * @tparam T template value.
   * @param rows number of rows.
   * @param cols number of columns
   * @return A identity matrix.
   */
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
  /**
   * @brief function to print a vector.
   * @tparam T template value.
   * @param vect vector to print.
   */
  template<typename T>
  static void vector_show(const std::vector<T> &vect){
    size_t size = vect.size();
    for(size_t i = 0; i< size ; ++i){
        std::cout << vect[i] << " " ;
      }
    std::cout << "\n";
  }
}

#endif