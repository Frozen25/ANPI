/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Crisptofer Fernandez
 * @Editor: Alexis Gavriel
 * @Date  : 07.10.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>
#include "Utilities.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"
#include <iostream>
#include "Intrinsics.hpp"
#include "IntrinsicsM.hpp"

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
    for (size_t col = 0; col < A.cols(); col++) {
      //max value of the current column
      T max = T(0);
      //index of the max value at current column
      size_t max_index = col;


      //In the following loop the max value and max index of column 'col' are setted
      for(size_t row = col; row < A.rows();row++){
        //searching the max absolute value
        if (std::abs(max) < std::abs(LU[row][col])){
          //setting the max value and max index
          max = LU[row][col];
          max_index = row;
        }
      }
      //We can swap the current row and the row which we have found the max value of the
      //current col. Remember, the current col is also the current row.
      
      pivot(LU,col,0,col,permut);  

      //Setting LU for each row, where row > col.
      for(size_t row = col+1; row < A.rows();row++){
        //verifing if the current element is zero.
        if (std::abs(LU[col][col]) == T(0))
          throw anpi::Exception("Division by zero detected, LU matrix couldn't be created");

        //Making the factor.
        T factor = LU[row][col]/LU[col][col];


        //U part of LU matrix
        for (size_t k = col+1; k < A.cols(); k++) {
          LU[row][k] -= LU[col][k]*factor;
        }
        //L part of U matrix
        LU[row][col] = factor;

      }


    }
        
      
    
  }
  

  namespace simd{
    #ifdef ANPI_ENABLE_SIMD
    #ifdef __AVX__


    template<typename T>
    void unpackDoolittleSIMD(const Matrix<T>& LU,
                         Matrix<T>& L,
                         Matrix<T>& U) {

        if (LU.rows() != LU.cols()) throw anpi::Exception("Matrix is not a square!");
        
        size_t n = LU.cols();     //Columns number;
        L.allocate(n, n);
        U.allocate(n, n);

        L.fill(T(0));
        U.fill(T(0));                              

        for (size_t i = 0; i < n; ++i) {             //Rows Iterator;

            for (size_t j = 0; j < n; ++j) {         //Columns Iterator;

                if(j>i) {                            //If current position is above the diagonal;
                    U[i][j] = LU[i][j];              //Add number from LU to U in the correspond space;                  
                }
                else if(j==i) {                      //If current position is in the diagonal;
                    U[i][j] = LU[i][j];              //Setting U's diagonal with LU's number;
                    L[i][j] = T(1);                     //Setting in one the L's diagonal;
                }
                else {                               //If current position is under the diagonal;                  
                    L[i][j] = LU[i][j];              //dd number from LU to U in the correspond space;
                }
            }

        }
      
    }///unpackDoolittleSIMD





    template<typename T,typename regType>
    void luDoolittleSIMD(const Matrix<T>& A,
                     Matrix<T>& LU,
                     std::vector<size_t>& permut) {

        if (A.rows() != A.cols()) throw anpi::Exception("Matrix is not a square!");
        
        size_t n = A.rows();                         //Rows number and Columns number, counter

        permut = std::vector<size_t> (n,0);    //Initialize the permutation vector
        for(size_t i = 0; i<permut.size(); ++i){      //Add elements to permutation vector
            permut[i] = i;
        }

        LU = A;
        //std::cout << "using luDoolittleSIMD" << std::endl;

        
        T factor = T(0);                                          //factor = factors for the lower triangular matrix

        const size_t tentries = A.rows() * A.dcols();     //Takes amount of entries in matrix A
        const size_t blocks =
                (tentries * sizeof(T) + (sizeof(regType) - 1)) / sizeof(regType); //Calculates amount of blocks
        const size_t blocks_in_row = blocks / A.rows();           //amount of blocks per row
        const size_t data_in_block = A.dcols() / blocks_in_row;   //amount of values per block

        
        std::vector<T> factor_vect (data_in_block, factor);   //fills vector with factor
       
        regType* factor_vector = reinterpret_cast<regType*>(&factor_vect);    //casts the pointer of the vector to register
        
        regType *start_LU = reinterpret_cast<regType*>(LU.data());    //pointer to the begining of the matrix
        //regType *end_LU   = start_LU + blocks;


        std::cout<< "starting SIMD LU...\n";

        pivot(LU,0,0,0,permut);                       //Pivoting
        for (size_t k = 0; k < n-1; ++k) {            //Elimination Iterator

          regType *block_actual = start_LU;                      //pointer to first block of this row
          regType *block_factor = start_LU + k * blocks_in_row;  
          regType *row_back =  block_actual;
          //regType *row_next = block_actual + blocks_in_row;        //pointer to first block of next row
          //regType *row_end = row_next + blocks_in_row;        //pointer to end of next block

          for (size_t i = k+1; (i < n) ; ++i) {       //Rows Iterator

            factor = LU[i][k]/LU[k][k];             //calculates new factor
            factor_vect = std::vector<T> (data_in_block, factor);     //updates vector with new factor
            regType *block_actual = start_LU + i * blocks_in_row;  

            for (size_t j = k; j < n; j+= data_in_block ) {      //Cols iterator

              //in case we have leftover data that wont fit in a reg
              if (j + data_in_block >= n){
                for( ; j < n; ++j){
                  LU[i][j] -= factor*LU[k][j];
                }
              }else{
                int h = (int) (k - kata_in_block * (j - 1));
                if (h >= 0) {
                  int missing = data_in_block - h; 
                  int m = k + 1;
                  int count = missing + m;
                  for (; m < count - 1; ++m) {
                      LU[i][m] -= factor * LU(k, m);
                    }//for j
                }else{
                //calculates the new values per block
                *block_actual = mm_sub<T, typename avx_traits<T>::reg_type>(
                    *block_actual, mm_mul<T, typename avx_traits<T>::reg_type>(*factor_vector,*block_factor)); 
                std::cout << "using simd operations...\n";
                // nonSIMD version:  LU[i][j] = LU[i][j] - factor*LU[k][j];
                }
              }

              ++block_actual;
                
            }///for cols iterator

            LU[i][k] = factor;                          //Filling the lower triangular matrix
            row_back  += blocks_in_row;
            block_actual == row_back; 
            //row_next  = row_back + blocks_in_row;

          }///for rows iterator


          pivot(LU,k,0,k,permut);                    //Pivoting
        }///for elimination iterator
      
      }///luDoolittleSIMD



      #endif
    #endif


  }//namespace simd


}//namespace anpi
  
#endif

