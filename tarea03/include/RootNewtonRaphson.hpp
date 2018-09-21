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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking by means of the
   * Newton-Raphson method
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial root guess
   * 
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {

    const int maxi = std::numeric_limits<T>::digits;           // Max number of iterations it's going to do
    T xr = xi;                                                 // Initializes the value of xr
    T fr  = funct(xi);                                         // Initializes the value of f(xr)

    T dfx = ( funct(xr + eps)- funct(xr - eps) )/(T(2)*eps);         //Initializes the value of derivative
    T h   = fr/dfx;                                            //Initializes the value of f(xr)/f'(xr)
    
    T ea = T();              // Approximate error

    for (int i = maxi; i > 0; --i) {
        
        xi = xr;
        xr = xi - h;      // New estimated root
        fr = funct(xr);   // New value of the root in the y axis
        ea = xr - xi;     // New Error


        if (std::abs(ea) < eps) { // Returns the value if the precision has been achieved
            return xr;
        }

        dfx = ( funct(xr+eps)- funct(xr - eps) )/(T(2)*eps); //New value of the derivative

        // Avoids division by zero
        if (std::abs(dfx) > std::numeric_limits<T>::epsilon()) {
            h = fr/dfx;
        }

        

        
    }


    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif
