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

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the Brent's method.
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
  T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps){
    // TODO: Put your code in here!
      if(xu<=xl){                                                          //EVALUATE IF THE VALUES OF X ARE INVERTED
          throw anpi::Exception("INVALID INTERVAL") ;
      }

      if(funct(xl)*funct(xu)>0){                                           //EVALUATE IF THE ENDPOINTS HAVE THE SAME SIGN
          throw anpi::Exception("BOTH POINTS HAVE THE SAME SIGN") ;
      }

      const int maxi= std::numeric_limits<T>::digits;                      //GET THE MAXIMUM NUMBER OF ITERATIONS
      T a = xl;                                                            //INITIALIZE THE A POINT
      T b = xu;                                                            //INITIALIZE THE B POINT
      T c = xu;                                                            //INITIALIZE THE C POINT
      T d,e,min1,min2;
      T fa= funct(a);
      T fb = funct(b);
      T fc, p,q,r,s,tol1,xm;

      if((fb>T(0) && fa >T(0)) || (fa<T(0) && fb < T(0))){
          return 0;
      }
      fc = fb;
      for (int i = 1; i<= maxi;i++){

          if ((fb > T(0) && fc > T(0)) || (fb < T(0) && fc < T(0))) {      //EVALUATE IF THERE ARE TEO POINTS WITH THE SAME SIGN
              c=a;                                                         //RENAME A,B,C AND ADJUST BOUNDING INTERVAL
              fc=fa;
              e=d=b-a;
          }

          if(std::fabs(fc) < std::fabs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
          }

          tol1=T(2)*eps*std::fabs(b)*T(0.5);                               //CONVERGENCE CHECK
          xm =T(0.5)*(c-b);

          if(std::abs(xm)<= tol1 || fb == T(0)){
            return b;

          }
          if(std::abs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)){
              s = fb/fa;                                                  //ATTEMPT INVERSE QUADRATIC INTERPOLATION
              if(a == c){
                  p= T(2)*xm*s;
                  q = T(1) -s;
              }
              else{
                  q = fa/fc;
                  r = fb/fc;
                  p=s*(T(2)*xm*q*(q-r)-(b-a)*(r-T(1)));
                  q = (q -T(0)) *(r -T(0))*(s-T(0));
              }
              if(p>T(0)){
                  q = -q;
              }
              
              p = std::abs(p);
              min1 = T(3) * xm * q - std::abs(tol1*q);
              min2 = std::abs(e*q);

              if(T(2)*p < (min1 < min2 ? min1 : min2)){

                  e=d;                                              //ACCEPT INTERPOLATION
                  d=p/q;
              }

              else{
                  d = xm;                                           //INTERPOLATION FAILED, USE BISECTION
                  e=d;
              }
          }
          else{                                                     //BOUND DECREASING TOO SLOWLY, USE BISECTION
              d = xm;
              e=d;
          }
          a = b;                                                    //MOVE LAST BEST GUESS TO A
          fa=fb;

          if(std::fabs(d) > tol1){                                  //EVALUATE NEW TRIAL ROOT
              b+=d;
          }
          else{
              b += ((xm) >= T(0) ? std::abs(tol1) : -std::abs(tol1));
          }
          fb = funct(b);
      }
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }
}
  
#endif

