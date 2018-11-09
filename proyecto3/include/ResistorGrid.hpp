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
  


  class ResistorGrid
  {
  private:
    /// Matrix of the current equation system
    Matrix<float> A_;
    /// Vector of the current equation system
    std::vector<float> b_;
    /// Vector of the current equation system
    std::vector<float> c_;


    bool isolatedTop, isolatedLeft, isolatedBottom, isolatedRight;
    
  public:
    

    ResistorGrid() {}
    
    

    /**
     * Construct the grid from the given file
     * @return true if successful or false otherwise
     */
    bool build(const std::string &filename) {

      

      return true;
    }

    
    
    float TopBar(float x){
      return 0;
    }
    
    float BottomBar(float x){
      return 0;
    }
    
    float LeftBar(float x){
      return 0;
    }
    
    float RightBar(float x){
      return 0;
    }

    size_t old_index(size_t yj){
      return (yj-1)/2 +1;
    }


    template<typename T>
    bool fill_borders(Matrix<T>& A){
      size_t rows = A.rows();
      size_t cols = A.cols();



      size_t Aj = 0;
      size_t Ai = 0;

      for( Aj = 1; Aj<(cols-1); ++Aj){
        A[0][Aj] = TopBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

      for( Aj = 1; Aj<(cols-1); ++Aj){
        A[rows-1][Aj] = BottomBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

      for( Ai = 1; Ai<(rows-1); ++Ai){
        A[Ai][0] = LeftBar( ((float)(Ai))/((float)(rows-2)) - ((float)(cols-2))/2);
      }

      for( Ai = 1; Ai<(rows-1); ++Aj){
        A[Ai][cols-1] = RightBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }
    
    }

    //this function creates a new matrix Y, which is double the size of the original inside its borders
    //  and fills it with the old matrix data, for every element in A there are 4 elements in Y
    template<typename T>
    void scale_matrix(const Matrix<T>&  A, Matrix<T>&  Y) {
      size_t rows = A.cols();
      size_t cols = A.rows();

      // creates a matrix with a size equal to the double of the amount of elements outside of the border
      Y.allocate(((rows-2)*2+2),((cols-2)*2+2));


      fill_borders(Y);

      size_t Ai = 0;
      size_t Aj = 0;
      size_t yi = 0;
      size_t yj = 0;
      double Aij_L = 0.0;
      double Aij_T = 0.0;
      double Aij_R = 0.0;
      double Aij_D = 0.0  ;

      for (yi = 1; yi < rows-1; ++yi){
        

        for( yj = 1; yj < cols-1; ++yj ){


        ////////////////////////////////////////////////////////////
        ///        I N T E R N A L   B L O C K 
        ////////////////////////////////////////////////////////////

          if( yi>1 && yi<(rows-1)  && yj>1 && yj<(cols-1) ){
          
            // maps the old index to the corresponding index in the new matrix
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
          
            //stores old values in new matrix
            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
        }


        ////////////////////////////////////////////////////////////
        ///             C O R N E R S 
        ////////////////////////////////////////////////////////////

        ///////////////////
        //  TOP LEFT
        ///////////////////
        else if ( (yi==1)&&(yj==1) ){
          if (isolatedTop && isolatedLeft ){
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = ( Aij_R + Aij_D)/2;
          }
          else if (isolatedTop){
            Aij_L = Y[1][0];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_L + Aij_R + Aij_D)/3;
          }
          else if (isolatedLeft){
            Aij_T = Y[0][1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_T + Aij_R + Aij_D)/3;
          }
          else{

            Aij_L = Y[1][0];
            Aij_T = Y[0][1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }

        ///////////////////
        //  TOP RIGHT
        ///////////////////
        else if ( (yi==1)&&(yj==(cols-2)) ){
          if (isolatedTop && isolatedRight ){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = ( Aij_L + Aij_D)/2;
          }
          else if (isolatedTop){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_R = Y[1][cols-1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_L + Aij_R + Aij_D)/3;
          }
          else if (isolatedRight){
            Aij_T = T[0][yj];
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_T + Aij_L + Aij_D)/3;
          }
          else{

            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = T[0][yj];
            Aij_R = Y[1][cols-1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }

        ///////////////////
        //  BOTTOM LEFT
        ///////////////////
        else if ( (yi==(rows-2))&&(yj==1) ){
          if (isolatedBottom && isolatedLeft ){
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Y[yi][yj] = ( Aij_T + Aij_R)/2;
          }
          else if (isolatedBottom){
            Aij_L = Y[yi][0];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Y[yi][yj] = (Aij_L + Aij_T + Aij_R)/3;
          }
          else if (isolatedLeft){
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = Y[1][rows-1];
            Y[yi][yj] = (Aij_T + Aij_R + Aij_D)/3;
          }
          else{

            Aij_L = Y[yi][0];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = Y[1][rows-1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }

        ///////////////////
        //  BOTTOM RIGHT
        ///////////////////
        else if ( (yi==(rows-2))&&(yj==(cols-2)) ){
          if (isolatedBottom && isolatedRight ){
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Y[yi][yj] = ( Aij_T + Aij_L)/2;
          }
          else if (isolatedBottom){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = Y[yi][cols-1];
            Y[yi][yj] = (Aij_L + Aij_T + Aij_R)/3;
          }
          else if (isolatedRight){
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_D = Y[yi+1][yj];
            Y[yi][yj] = (Aij_T + Aij_L + Aij_D)/3;
          }
          else{

            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = Y[yi][cols-1];
            Aij_D = Y[yi+1][yj];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }


        ////////////////////////////////////////////////////////////
        ///             B O R D E R S
        ////////////////////////////////////////////////////////////


        ///////////////////
        //  TOP BORDER
        ///////////////////
        else if ( yi==1 ){
          if (isolatedTop){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_L + Aij_R + Aij_D)/3;          
          }
          else{

            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = Y[0][yi];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }
        ///////////////////
        //  BOTTOM BORDER
        ///////////////////
        else if ( yi==(rows-2) ){
          if (isolatedBottom){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_L + Aij_R + Aij_T)/3;          
          }
          else{

            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = Y[yi+1][yj];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }
        ///////////////////
        //  LEFT BORDER
        ///////////////////
        else if ( yj==1 ){
          if (isolatedLeft){
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_T + Aij_R + Aij_D)/3;          
          }
          else{

            Aij_L = Y[yi][0];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }
        ///////////////////
        //  RIGHT BORDER
        ///////////////////
        else if ( yj==(cols-2) ){
          if (isolatedLeft){
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_T + Aij_L + Aij_D)/3;          
          }
          else{

            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = Y[yi][yj+1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }


        

        }/// for yj
      } ///  for yi
    }   /// scale matrix

    

    
  };
}


#endif //PROYECTO2_RESISTORGRID_HPP
