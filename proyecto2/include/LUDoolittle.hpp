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



      /*
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
        const size_t cols_in_block = A.dcols() / blocks_in_row;   //amount of values per block

        
        //std::vector<T,anpi::aligned_allocator<T> > factor_vector (cols_in_block, factor);   //fills vector with factor
        Matrix<T> factor_vector;//This matrix is used to create pointer to a vector a lot of times
        factor_vector.allocate(1, cols_in_block);
        regType *factor_reg = reinterpret_cast<regType*>(factor_vector.data());  //casts the pointer of the vector to register
        
        regType *start_LU = reinterpret_cast<regType*>(LU.data());    //pointer to the begining of the matrix
        //regType *end_LU   = start_LU + blocks;


        pivot(LU,0,0,0,permut);                       //Pivoting
        for (size_t col = 0; col < n-1; ++col) {            //Elimination Iterator

          regType *block_col_actual = reinterpret_cast<regType*>(&LU[col][col]);                     //pointer to first block of this row
          regType *block_col_backup = block_col_actual;
          regType *block_row_actual;
          
          //regType *row_next = block_col_actual + blocks_in_row;        //pointer to first block of next row
          //regType *row_end = row_next + blocks_in_row;        //pointer to end of next block

          for (size_t row = col+1; row < n ; ++row) {       //Rows Iterator

            block_row_actual = reinterpret_cast<regType*>(&LU[row][row]);
            
            
            factor = LU[row][col]/LU[col][col];             //calculates new factor
            //factor_vector = std::vector<T,anpi::aligned_allocator<T> > (cols_in_block, factor);     //updates vector with new factor
            factor_vector.fill(factor);
            factor_reg = reinterpret_cast<regType*>(factor_vector.data());

            regType data1;
            //regType data2;
            for (size_t j = col; j < LU.cols(); j+= cols_in_block ) {      //Cols iterator

                //calculates the new values per block
                
                data1 = mm_mul<T, regType>(*factor_reg,*block_col_actual);
                *block_row_actual =  mm_sub<T, regType>(*block_row_actual, data1);
                  *block_row_actual = data1;

                  *block_row_actual = mm_sub<T, regType>(
                    *block_row_actual, mm_mul<T, regType>(*factor_reg,*block_col_actual)); 
                
                // nonSIMD version:  LU[row][j] = LU[row][j] - factor*LU[k][j];
              ++block_row_actual;  
              ++block_col_actual;
                
            }///for cols iterator

            LU[row][col] = factor;                          //Filling the lower triangular matrix
            block_col_actual = block_col_backup;
            
            //row_next  = row_back + blocks_in_row;

          }///for rows iterator


          pivot(LU,col,0,col,permut);                    //Pivoting
        }///for elimination iterator
      
      }///luDoolittleSIMD
      */


    // ludolittle for simd
    template<typename T,typename regType>
    inline void luDoolittleSIMD(const Matrix<T>& A,
                                 Matrix<T>& LU,
                               std::vector<size_t>& permut) {
      if(A.rows() != A.cols()) throw anpi::Exception("Matrix is not square"); 
      
      LU = Matrix<T>(A); 

      const size_t tentries = A.rows()*A.dcols();         // amount of values in matrix
      const size_t blocks = (tentries*sizeof(T) + (sizeof(regType)-1) )/sizeof(regType);      //amount of blocks in matrix 


      const size_t blocks_in_row = blocks/A.rows();       //amount off blocks in a row
      const size_t cols_in_block = A.dcols()/blocks_in_row;       ///amount of values in a block

      Matrix<T> factor_vector;                        //a matrix is used to store the factor as many times as values in a block
      factor_vector.allocate(1, cols_in_block);       //this is due to alignment of data
      regType* factor_reg  = reinterpret_cast<regType*>(factor_vector.data());   //casts to a register

      permut = std::vector<size_t>(A.rows(), 0); 
      for(size_t x = 0; x < permut.size(); ++x){ 
          permut[x] = x;
      }


      regType* LU_start  = reinterpret_cast<regType*>(LU.data());     //saves the begining of the matrix
      regType* LU_end  = LU_start + blocks;                           //saves the end of the matrix


      for(size_t col = 0; col < LU.rows() - 1; ++col){                      //calculates factos of column col
          pivot(LU, col, 0, col, permut);                                   ///pivoting
          regType* block_col_actual = LU_start + col * blocks_in_row;       ///saves the actual row
          regType* actual_col_backup = block_col_actual;                    //makes a backup of the row
          regType* block_row_actual = block_col_actual + blocks_in_row;     //saves the actual row to iterate
          regType* end_next_row = block_row_actual + blocks_in_row;         //saves the end of the row to iterate

          T factor;
          for(size_t row = col + 1; block_row_actual != LU_end;++row){       //calculates factors of column col and row row

              if(std::abs(LU(col, col)) < std::numeric_limits<T>::epsilon()){ 
                  throw anpi::Exception("Division by zero! Matrix cannot be decomposed by luDoolittle");
              }
              factor = LU(row, col)/LU(col, col);         //factor that multiplies every value
              LU[row][col] = factor;                      //stores the factor
              
              for(size_t k = 1;block_row_actual!=end_next_row;++block_row_actual){ 
                                                      //iterates on the row, updates it by the formula:
                                                                  // row = row - factor*row2
                                                                  // where row2 is the LU[col] row
                  
                  if(cols_in_block * k - 1 >= col){       //checks if the values overlap into the U matrix
                      factor_vector.fill(factor); 
                      
                      
                      int over = (int)(col - cols_in_block * (k - 1)); ///checks to prevent destroying past factors
                      for(int count = 0; count <= over; ++count){ 
                          factor_vector[0][count] = T(0);             ///to prevent it, fills the factor vector with 0
                      }
                      ///NON SIMD: LU[row][k] -= LU[col][k]*factor;                        
                      *block_row_actual = mm_sub<T>(*block_row_actual, mm_mul<T,regType>(*factor_reg, *block_col_actual)); 
                                                                                                                
                  }
                  ++k;
                  ++block_col_actual;
                  
              }
              
              end_next_row += blocks_in_row;
              block_col_actual = actual_col_backup;
              
          }///END row iterarion
        }///END col iteration
      }///END DoolittleSimd
  

  



      #endif
    #endif


  }//namespace simd


}//namespace anpi
  
#endif

