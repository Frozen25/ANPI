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
#include "omp.h"

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

    /// Epsilon used to calculate convergence
    double epsilon = 1.0001f;

    // Convergence flag
    bool convergenceFlag = false;

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
      return 8.0f;
    }


    /**
     * @brief [calculates the corresponding matching index on the new matrix (Y)]
     * @details [this is the mapping of the index from the small matrix to the bigger matrix]
     * 
     * @param yj [index in the matrix A]
     * @return [index in the matrix Y]
     */
    size_t new_index(size_t Ax){ //TODO: check method
      return (Ax-1)*2 +1;
    }

    /**
     * @brief [calculates the corresponding matching index on the old matrix (A)]
     * @details [this is the inverse of mapping the index from the big matrix to the smaller matrix]
     * 
     * @param yj [index in the matrix Y]
     * @return [index in the matrix A]
     */
    size_t old_index(size_t Yx){ //TODO: check method
      return (Yx-1)/2 +1;
    }

    /**
     * @brief fills the borders of the matrix with the corresponding values
     * @details [long description]
     * 
     * @param A A matrix that contains the thermal information.
     * @tparam T template value.
     */
    template<typename T>
    void fill_borders(Matrix<T>& Mat  ){
      size_t rows = Mat.rows();
      size_t cols = Mat.cols();      
      size_t Matj = 0;
      size_t Mati = 0;

      /// The middle point of the block being filled is normalized with: (float)(2*VALUE-1)/(2*(cols-2))
          
      /// This fills the Top Bar
      for( Matj = 1; Matj<(cols-1); ++Matj){
        Mat[0][Matj] = TopBar( (float)(2*Matj-1)/(2*(cols-2)), rows );
      }

      /// This fills the Bottom Bar
      for( Matj = 1; Matj<(cols-1); ++Matj){
        Mat[rows-1][Matj] = BottomBar( (float)(2*Matj-1)/(2*(cols-2)) );
      }

      /// This fills the Left Bar
      for( Mati = 1; Mati<(rows-1); ++Mati){
        Mat[Mati][0] = LeftBar( (float)(2*Matj-1)/(2*(cols-2)) );
      }

      /// This fills the Right Bar
      for( Mati = 1; Mati<(rows-1); ++Mati){
        Mat[Mati][cols-1] = RightBar( (float)(2*Matj-1)/(2*(cols-2)) );
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

      
      /// Default Initialized values
      double Aij_L = 0.0f;
      double Aij_T = 0.0f;
      double Aij_R = 0.0f; 
      double Aij_D = 0.0f; 

      bool convergenceFlagParallel = true;
      
      #pragma omp parallel for default(none) private(Aij_L, Aij_T, Aij_R, Aij_D) shared(A , Y, yrows, ycols, convergenceFlagParallel)
      for (size_t yi = 1; yi < yrows-1; ++yi){
        for(size_t yj = 1; yj < ycols-1; ++yj ){
          
          
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
          else if ( (yi==(yrows-2))&&(yj==(ycols-2)) ){
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
          else if ( yj==(ycols-2) ){
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

          /// checks convergence, if it converges, the final result is true
          if(std::abs ( Y[yi][yj] - A[(yi-1)/2+1][(yj-1)/2+1] ) > epsilon  ){
            convergenceFlagParallel = false;
          }

        }// for yj
      }// for yi
      convergenceFlag = convergenceFlagParallel;
    }// scale matrix



    /// Default value of max iterations is set to 15.
    template<typename T>
    void calculatePlate(Matrix<T>&  A, Matrix<T>&  Y , double eps, size_t maxIterations = 10){
      
      if (eps)
        epsilon = eps;

      //anpi::Matrix<double> A;
      //anpi::Matrix<double> Y;
      A.allocate(3,3);
      fill_borders(A);
      A[1][1] = ( A[0][1] + A[1][0] + A[1][2] + A[2][1] )/4;
      Y = A;
      for (size_t k = 0; (k < maxIterations && (!convergenceFlag)); ++k){
        A = Y;
        scale_matrix(A,Y);        
        std::cout << "iteration: "<< k << "\tnumber of rows = "<< A.rows() << std::endl;          
      }


      /// changes the output of printing a boolean "1" to "true", and a "0" to "false"
      std::cout << std::boolalpha;
      std::cout << "convergence: " <<convergenceFlag << std::endl;

      /// Function used to save the matrix in a file called  matrix.txt
      //matrix_show_file(Y);

    }






  };
}

#endif //PROYECTO3_THERMALPLATE_HPP
