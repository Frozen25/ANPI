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

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the bisection method.
   *
   * @param funct a std::function of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    if (xu < xl) { // Throws exception if the intervals are inverted

      throw anpi::Exception("Inverted intervals");

    } else {

        T fl = funct(xl);        // Value of the lower limit in the y axis. Saves calculations
        T fu = funct(xu);        // Value of the upper limit in the y axis. Saves calculations

        if (fu*fl > 0) { // Throws exception if funct(xu) and funct(xl) have the same sign

            throw anpi::Exception("Unenclosed root");

        } else { // Solves the problem

            const int maxi = std::numeric_limits<T>::digits; // Max number of iterations it's going to do

            T xr = xl;                // Initializes the value with something valid
            T ea = T();              // Approximate error

            for (int i = maxi; i > 0; --i) {
                T xrold(xr);       // Old root, necessary to calculate the error
                xr = (xl + xu) / T(2);  // New estimated root
                T fr = funct(xr);   // value of the root in the y axis

                // Avoids division by zero
                if (std::abs(xr) > std::numeric_limits<T>::epsilon()) {
                    ea = std::abs((xr - xrold) / xr);
                }

                T cond = fl * fr;  // value to know who will be the new end
                                            // depending on the sign of the y axis

                if (cond < T(0)) {
                    xu = xr;              // Different sign the root and lower, so changes the upper limit
                } else if (cond > T(0)) {
                    xl = xr;              // Same sign the root and lower, so changes the lower limit
                    fl = fr;         // updates the value in the y axis because it's pre-image changed
                } else {
                    ea = T(0);           // No error! => some value is zero
                    xr = (std::abs(fl) < std::numeric_limits<T>::epsilon())
                             ? xl : xr;   // f_lower == 0
                }

                if (ea < eps) { // Returns the value if the precision has been achieved
                    return xr;
                }
            }
        }
    }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

