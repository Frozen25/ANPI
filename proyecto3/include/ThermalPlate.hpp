//
// Created by ubuntu on 01/10/18.
//

#ifndef PROYECTO3_THERMALPLATE_HPP
#define PROYECTO3_THERMALPLATE_HPP

#include <cstdlib>
#include <iostream>
#include <string>

#include <AnpiConfig.hpp>
#include <Spline.h>

#include <Matrix.hpp>
#include <Exception.hpp>
#include "Solver.hpp"
#include "Utilities.hpp"


namespace anpi{

/** class ThermalPlate
 *  @brief thermal plate class
 *  @details This class is an abstraction of a thermal plate,
 *  is used to make the temperature gradient using the liebman
 *  method.
 */
  class ThermalPlate
  {
  private:
    /// Matrix of the current equation system
    Matrix<float> A_;
    /// Vector of the current equation system
    std::vector<float> b_;
    /// Vector of the current equation system
    std::vector<float> c_;

    /// Flags that indicates if a side of the plate is isolated
    bool isolatedTop = false, 
         isolatedLeft = false, 
         isolatedBottom = false, 
         isolatedRight = false;
    
    /// Vectors with the temperatures of the sides
    std::vector<double> top, bottom, left, right;
    
    /// Flag to show the thermal flow of the plate
    bool showThermalFlow{};
    
    /// Flag to show no visuals (overrides every other visual parameter)
    bool noVisuals{};
    
    /// Size of the desired grid to show (it's square)
    int grid{}; // Each "grid" pixels, an arrow of the flow is shown
    int horizontal{};
    int vertical{};
    
  public:
    
    /// Default constructor
    ThermalPlate() = default;

  /**
   * Constructor of the class.
   * @tparam T type of the matrix
   * @param isolatedTop used to determine if the top of the matrix is isolated.
   * @param isolatedBottom used to determine if the bottom of the matrix is isolated.
   * @param isolatedLeft used to determine if the left of the matrix is isolated.
   * @param isolatedRight used to determine if the right of the matrix is isolated.
   * @param top This vector describe the thermal profile in the top of the thermal plate.
   * @param bottom This vector describe the thermal profile in the bottom of the thermal plate.
   * @param left This vector describe the thermal profile in the left of the thermal plate.
   * @param right This vector describe the thermal profile in the right of the thermal plate.
   * @param horizontal horizontal value.
   * @param vertical vertical value.
   * @param grid grid used to show the thermal flow.
   * @param noVisuals Flag to show no visuals.
   * @param showThermalFlow flag used to visualize the thermal flow.
   */
    ThermalPlate(bool isolatedTop, bool isolatedBottom, bool isolatedLeft, bool isolatedRight,
        std::vector<double> &top, std::vector<double> &bottom, std::vector<double> &left, std::vector<double> &right,
        bool showThermalFlow, bool noVisuals, int grid, int horizontal, int vertical) {
      
      this->isolatedTop = isolatedTop;
      this->isolatedBottom = isolatedBottom;
      this->isolatedLeft = isolatedLeft;
      this->isolatedRight = isolatedRight;
      
      this->top = top;
      this->bottom = bottom;
      this->left = left;
      this->right = right;
      
      this->showThermalFlow = showThermalFlow;
      this->noVisuals = noVisuals;
      this->grid = grid;
      this->horizontal = horizontal;
      this->vertical = vertical;
    }
    
    double cubicSpline(std::vector<double> *border, double &position){
      std::vector<double> x; // The normalized x position of each temperature in the vector border
      size_t numTemperature = border->size(); // Size of the vector border
      
      for(size_t i = 0; i < numTemperature; ++i){
        x.push_back(static_cast<double>(i)/(numTemperature-1));
      }
  
      tk::spline s;
      s.set_points(x, *border);    // currently it is required that X is already sorted
  
      return s(position);
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

    size_t old_index(size_t yj){ //TODO: check method
      return (yj-1)/2 +1;
    }

    /**
     * @brief Used to fill the borders with the thermal profile that user gives.
     * @details The method uses a matrix to fill its edges with a thermal profile,
     * then it is used to expand this temperature to the rest of the matrix.
     * 
     * @param A A matrix that contains the thermal information.
     * @tparam T template value.
     */
    template<typename T>
    void fill_borders(Matrix<T>& A){
      size_t rows = A.rows();
      size_t cols = A.cols();
      
      size_t Aj = 0;
      size_t Ai = 0;

      for( Aj = 1; Aj<(cols); ++Aj){
        A[0][Aj] = TopBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

      for( Aj = 1; Aj<(cols); ++Aj){
        A[rows-1][Aj] = BottomBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

      for( Ai = 1; Ai<(rows); ++Ai){
        A[Ai][0] = LeftBar( ((float)(Ai))/((float)(rows-2)) - ((float)(cols-2))/2);
      }

      for( Ai = 1; Ai<(rows); ++Ai){
        A[Ai][cols-1] = RightBar( ((float)(Aj))/((float)(cols-2)) - ((float)(cols-2))/2);
      }

    }

    /**
     * Creates a new matrix Y, which is double the size of the original inside its borders
     * and fills it with the old matrix data, for every element in A there are 4 elements in Y
     * @tparam T type of the matrix
     * @param A old matrix (small one)
     * @param Y new matrix (big one)
     */
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
      double Aij_R = 0.0; //TODO: check value
      double Aij_D = 0.0; //TODO: check value

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
          ///  TOP LEFT
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
          ///  TOP RIGHT
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
          ///  BOTTOM LEFT
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
          ///  BOTTOM RIGHT
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
          ///  TOP BORDER
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
          ///  BOTTOM BORDER
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
          ///  LEFT BORDER
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
          ///  RIGHT BORDER
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
        }// for yj
      }// for yi
    }// scale matrix
  };
}

#endif //PROYECTO3_THERMALPLATE_HPP
