//
// Created by ubuntu on 01/10/18.
//

#ifndef PROYECTO2_THERMALPLATE_HPP
#define PROYECTO2_THERMALPLATE_HPP

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
  


  class ThermalPlate
  {
  private:
    /// Matrix of the current equation system
    Matrix<float> A_;
    /// Vector of the current equation system
    std::vector<float> b_;
    /// Vector of the current equation system
    std::vector<float> c_;


    bool isolatedTop = false, 
         isolatedLeft = false, 
         isolatedBottom = false, 
         isolatedRight = false;
    
  public:
    

    ThermalPlate() {}
    
    

    /**
     * Construct the grid from the given file
     * @return true if successful or false otherwise
     */
    bool build(const std::string &filename) {

      

      return true;
    }

    
    
    float TopBar(float x){
      return 3.0f;
    }
    
    float BottomBar(float x){
      return 2.0f;
    }
    
    float LeftBar(float x){
      return -1.0f;
    }
    
    float RightBar(float x){
      return 4.0f;
    }

    size_t old_index(size_t yj){
      return (yj-1)/2 +1;
    }

    /**
     * @brief [brief description]
     * @details [long description]
     * 
     * @param A [description]
     * @tparam T [description]
     */
    template<typename T>
    void fill_borders(Matrix<T>& A  ){
      size_t rows = A.rows();
      size_t cols = A.cols();



      size_t Aj = 0;
      size_t Ai = 0;

      //printf("Filling TOP BAR\n");
      
      for( Aj = 1; Aj<(cols); ++Aj){
        A[0][Aj] = TopBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }


      //printf("Filling BOTTOM BAR\n");
      for( Aj = 1; Aj<(cols); ++Aj){
        A[rows-1][Aj] = BottomBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

      //printf("Filling LEFT BAR\n");
      for( Ai = 1; Ai<(rows); ++Ai){
        A[Ai][0] = LeftBar( ((float)(Ai))/((float)(rows-2)) - ((float)(cols-2))/2);
      }

      //printf("Filling RIGHT BAR\n");
      for( Ai = 1; Ai<(rows); ++Ai){
        A[Ai][cols-1] = RightBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

      //printf("Filled all bars\n");

    
    }

    //this function creates a new matrix Y, which is double the size of the original inside its borders
    //  and fills it with the old matrix data, for every element in A there are 4 elements in Y
    template<typename T>
    void scale_matrix( Matrix<T>&  A, Matrix<T>&  Y) {
      size_t rows = A.cols();
      size_t cols = A.rows();

      // creates a matrix with a size equal to the double of the amount of elements outside of the border
      Y.allocate(((rows-2)*2+2),((cols-2)*2+2));

      size_t yrows = (rows-2)*2+2;
      size_t ycols = (cols-2)*2+2;

      fill_borders(Y);

      size_t Ai = 0;
      size_t Aj = 0;
      size_t yi = 0;
      size_t yj = 0;
      double Aij_L = 0.0;
      double Aij_T = 0.0;
      double Aij_R = 0.0;
      double Aij_D = 0.0;
      
      

      for (yi = 1; yi < yrows-1; ++yi){
        

        for( yj = 1; yj < ycols-1; ++yj ){


        ////////////////////////////////////////////////////////////
        ///        I N T E R N A L   B L O C K 
        ////////////////////////////////////////////////////////////

          if( yi>1 && yi<(yrows-1)  && yj>1 && yj<(ycols-1) ){
          
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
          printf("%f , %f, %f, %f\n",Aij_L, Aij_T, Aij_R, Aij_D );
        }

        ///////////////////
        //  TOP RIGHT
        ///////////////////
        else if ( (yi==1)&&(yj==(ycols-2)) ){
          if (isolatedTop && isolatedRight ){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = ( Aij_L + Aij_D)/2;
          }
          else if (isolatedTop){
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_R = Y[1][yj+1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_L + Aij_R + Aij_D)/3;
          }
          else if (isolatedRight){
            Aij_T = Y[0][yj];
            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];
            Y[yi][yj] = (Aij_T + Aij_L + Aij_D)/3;
          }
          else{

            Aij_L = A[(yi-1)/2 +1][(yj-2)/2 +1];
            Aij_T = Y[0][yj];
            Aij_R = Y[1][yj+1];
            Aij_D = A[(yi  )/2 +1][(yj-1)/2 +1];

            Y[yi][yj] = (Aij_L + Aij_T + Aij_R + Aij_D)/4;
          }
        }

        ///////////////////
        //  BOTTOM LEFT
        ///////////////////
        else if ( (yi==(yrows-1))&&(yj==1) ){
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
            Aij_D = Y[yi+1][1];
            Y[yi][yj] = (Aij_T + Aij_R + Aij_D)/3;
          }
          else{

            Aij_L = Y[yi][0];
            Aij_T = A[(yi-2)/2 +1][(yj-1)/2 +1];
            Aij_R = A[(yi-1)/2 +1][(yj  )/2 +1];
            Aij_D = Y[yi+1][1];

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
            Aij_R = Y[yi][yj+1];
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
            Aij_R = Y[yi][yj+1];
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
        else if ( yi==(yrows-2) ){
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


#endif //PROYECTO2_THERMALPLATE_HPP
