/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

namespace anpi {
  
  /**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
  template<typename T>
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T eps) {

      const int maxi = std::numeric_limits<T>::digits;           // Max number of iterations it's going to do

      T xp;
      T xc;

      if (std::abs(funct(xi)) < std::abs(funct(xii))){           //Choose the lowest function value as the current x
          xp = xii;                                            //Initialize the previous value of x
          xc = xi;                                              //Initialize the current value fo x
      } else {
          xp = xi;                                             //Initialize the previous value of x
          xc = xii;                                            //Initialize the current value fo x
      }

      T ea;                                                      //Initialize the error value
      T dfx;                                                     //Initialize f(x-1) - f(x)
      T temp;                                                    //Initialize a temporary variable

      for(int i = maxi; i > 0; --i){

          ea = xp - xc;                                          //New error

          if (std::abs(ea) < eps) {                              // Returns the value if the precision has been achieved
              return xc;
          }
          dfx = funct(xp) - funct(xc);

          // Avoids division by zero
          if (std::abs(dfx) > std::numeric_limits<T>::epsilon()) {
              temp = xc;                                         //temporary storage to current value of x
              xc = xc - (funct(xc)*ea)/dfx;                      //Calculate the new value for the current x
              xp = temp;                                         //Assign the value of current x to previous x
          }
      }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

