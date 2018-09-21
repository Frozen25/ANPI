/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 04.08.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_RIDDER_HPP
#define ANPI_ROOT_RIDDER_HPP

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
  T rootRidder(const std::function<T(T)>& funct,T xi,T xii,const T eps) {
    if (xii < xi) { // Throws exception if the intervals are inverted
          throw anpi::Exception("Inverted intervals");
    }
    else{
      T fi = funct(xi);     //f(xlower)
      T fii = funct(xii);   //f(xupper)

      const int maxi = std::numeric_limits<T>::digits;           // Max number of iterations it's going to do

      if (fi*fii < T(0)){   //verifies if f(xlow) and f(xupper) have different sign

        T xm, fm, fr, s, ea;                  // xm = xmid , f(xmid) , s = geometric median of middle point and upper and lower ends, ea = Approximate Relative Error.
        T xr = T(0);                                               //  Initializes the value with something valid

        for (int i = maxi; i > 0; --i) {      //iterates for the number of max iterations

          T xrold(xr);                                          // Old root, necessary to calculate the error

          xm = 0.5 * (xi + xii);
          fm = funct(xm);         //evaluates function in the middle point of xi and xii
          s = std::sqrt(fm * fm - fi * fii);    //calculates a geometric average of the function in the middle point and the upper and lower ends

          xr = xm + (xm - xi) * ((fi >= fii ? T(1) : T(-1)) * fm / s);  //new possible value
          fr = funct(xr);   //evaluates function in new possible value

          // Avoids division by zero
          if (std::abs(xr) > std::numeric_limits<T>::epsilon()) {
            ea = std::abs((xr - xrold) / xr);   //new approximate error
          }
          //reassigns variables for next iteration
          if (fm * fr < T(0)) {
            xi = xm;
            fi = fm;
            xii = xr;
            fii = fr;
          } else if (fi * fr < T(0)) {
            xii = xr;
            fii = fr;
          } else if (fii * fr < T(0)) {
            xi = xr;
            fi = fr;
          }

          if (ea < eps) { // Returns the value if the precision has been achieved
            return xr;
          }
        }
        //if the root is not enclosed in the range
      } else{
        if ( abs (fi) == T(0)) {
          return xi;
        } else if ( abs (fii) == T(0)) {
          return xii;
        } else { // Throws exception if funct(xu) and funct(xl) have the same sign
          throw anpi::Exception("Unenclosed root");
        }
      }
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

