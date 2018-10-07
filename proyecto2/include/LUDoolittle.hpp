/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include "Utilities.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"
#include <iostream>

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {


  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of L
   * is composed by 1's.
   */
  template<typename T>
  void unpackDoolittle(const Matrix<T>& LU,
                       Matrix<T>& L,
                       Matrix<T>& U) {

      L = LU;
      U = LU;
      size_t n = LU.cols();                        //Columns number;

      if (LU.rows() != LU.cols()) throw anpi::Exception("Matrix is not a square!");

      for (size_t i = 0; i < n; ++i) {             //Rows Iterator;

          for (size_t j = 0; j < n; ++j) {         //Columns Iterator;

              if(j>i) {                            //If current position is above the diagonal;
                  U[i][j] = LU[i][j];              //Add number from LU to U in the correspond space;
                  L(i,j) = T(0);                      //Setting zero the space in the L matrix;
              }
              else if(j==i) {                      //If current position is in the diagonal;
                  U[i][j] = LU[i][j];              //Setting U's diagonal with LU's number;
                  L[i][j] = T(1);                     //Setting in one the L's diagonal;
              }
              else {                               //If current position is under the diagonal;
                  U[i][j] = T(0);                     //Setting zero the space in the U matrix;
                  L[i][j] = LU[i][j];              //dd number from LU to U in the correspond space;
              }
          }

      }
    //throw anpi::Exception("To be implemented yet");
  }
  
  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU. 
   *
   * The L matrix will have in the Doolittle's LU decomposition a
   * diagonal of 1'sc
   *
   * @param[in] A a square matrix 
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   */
  
  template<typename T,typename regType>
  void luDoolittleSIMD(const Matrix<T>& A,
                   Matrix<T>& LU,
                   std::vector<size_t>& permut) {

      if (A.rows() != A.cols()) throw anpi::Exception("Matrix is not a square!");
      permut = std::vector<size_t> (A.rows(),0);    //Initialize the permutation vector
      for(size_t i = 0; i<permut.size(); ++i){      //Add elements to permutation vector
          permut[i] = i;
      }

      LU = A;
      size_t n = LU.rows();                         //Rows number and Columns number, counter
      T f;                                          //f = factors for the lower triangular matrix

      regType element;
      
      pivot(LU,0,0,0,permut);                       //Pivoting
      for (size_t k = 0; k < n-1; ++k) {            //Column Iterator

          for (size_t i = k+1; i < n ; ++i) {       //Rows Iterator
              f = LU[i][k]/LU[k][k];

              for (size_t j = k; j < n; ++j) {      //Rows iterator for the elemental operation

                  LU[i][j] = LU[i][j] - f*LU[k][j];
              }
              LU[i][k] = f;                          //Filling the lower triangular matrix
          }
          pivot(LU,k,0,k,permut);                    //Pivoting
      }
    //throw anpi::Exception("To be implemented yet");
  }





  template<typename T>
  void luDoolittle(const Matrix<T>& A,
                   Matrix<T>& LU,
                   std::vector<size_t>& permut) {

      if (A.rows() != A.cols()) throw anpi::Exception("Matrix is not a square!");
      permut = std::vector<size_t> (A.rows(),0);    //Initialize the permutation vector
      for(size_t i = 0; i<permut.size(); ++i){      //Add elements to permutation vector
          permut[i] = i;
      }

      LU = A;
      size_t n = LU.rows();                         //Rows number and Columns number, counter
      T f;                                          //f = factors for the lower triangular matrix

      pivot(LU,0,0,0,permut);                       //Pivoting
      for (size_t k = 0; k < n-1; ++k) {            //Column Iterator

          for (size_t i = k+1; i < n ; ++i) {       //Rows Iterator
              f = LU[i][k]/LU[k][k];

              for (size_t j = k; j < n; ++j) {      //Rows iterator for the elemental operation

                  LU[i][j] = LU[i][j] - f*LU[k][j];
              }
              LU[i][k] = f;                          //Filling the lower triangular matrix
          }
          pivot(LU,k,0,k,permut);                    //Pivoting
      }
    //throw anpi::Exception("To be implemented yet");
  }

}
  
#endif

