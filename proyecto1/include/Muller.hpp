/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_MULLER_HPP
#define ANPI_MULLER_HPP

#include <vector>
#include <type_traits>

#include <PolynomialFormulaFormat.hpp>
#include <Deflation.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>


namespace anpi {

  /// Enumerator makes explicit polishing roots
  enum PolishEnum {
    DoNotPolish,
    PolishRoots
  };

  namespace bmt=boost::math::tools; // for polynomial

  /**
   * Finds a root of a given polynomial expression using muller's method
   * @param[in] poly polynomial to be analyzed for roots
   * @param[in] x point with number where to start to search for a root
   * @param[out] root stores the root found
   * @param[in] eps wanted error for the result
   *
   * @return root of the polynomial expression
   */
  template<class T,class U>
  typename std::enable_if<std::is_floating_point<U>::value,void>::type
  mullerRoots(const bmt::polynomial<T>& poly, const bmt::polynomial<T>& originalPoly,
              const U& x, U& root, const PolishEnum polish = DoNotPolish) {

    const int maxi = std::numeric_limits<U>::digits; // Max number of iterations it's going to do
    const U& eps = std::numeric_limits<U>::epsilon();

    //size of the polynomial is 0
    if (poly.degree() < 1){
      return;
    }

    bool found = false;                             // Bool to know if the root was found
    U unpolishRoot;

    U x0 = x;
    U x1 = x0 + U(1);
    U x2 = x1 + U(1);

    U f0 = bmt::evaluate_polynomial(&poly.data()[0], x0, poly.size());
    U f1 = bmt::evaluate_polynomial(&poly.data()[0], x1, poly.size());

    if (std::abs(f0) <= eps){                                              // Evaluates if the first and second numbers are the root
      unpolishRoot = x0;
      found = true;
    } else if(std::abs(f1) <= eps) {
      unpolishRoot = x1;
      found = true;
    }

    for(int i = 0; i < maxi; ++i){
      if(!(found)){
        U f2 = bmt::evaluate_polynomial(&poly.data()[0], x2, poly.size());

        if(std::abs(f2) <= eps){
          unpolishRoot = x2;
          found = true;
          break;
        }
        // Variables to resolve the muller's method
        U h0 = x1 - x0;
        U h1 = x2 - x1;

        U d0 = (f1 - f0) / h0;
        U d1 = (f2 - f1) / h1;

        U a = (d1 - d0) / (h1 - h0);
        U b = a * h1 + d1;
        U c = f2;

        // Because the method produces two roots, chooses the one with the same sign as b.
        U x3 = (b >= U(0)) ? x2 - ((U(-2) * c) / (b + sqrt(b * b - U(4) * a * c)))
                           : x2 - ((U(-2) * c) / (b - sqrt(b * b - U(4) * a * c)));

        x2 = x3;
        x1 = x2;
        x0 = x1;

        f1 = f2;
        f0 = f1;
      } else {
        break;
      }
    }
    //Polishes if necessary
    if(found && (polish == PolishRoots)){
      anpi::mullerRoots(originalPoly, originalPoly, unpolishRoot, root);
    } else if (found){
      root = unpolishRoot;
    }
  }

  template<class T,class U>
  typename std::enable_if<boost::is_complex<U>::value,void>::type
  mullerRoots(const bmt::polynomial<T>& poly, const bmt::polynomial<T>& originalPoly,
              const U& x, U& root, const PolishEnum polish = DoNotPolish) {

    //typedef to simplify instantiations of complex numbers
    typedef typename anpi::detail::inner_type<U>::type Utype;

    const int maxi = std::numeric_limits<Utype>::digits;       // Max number of iterations it's going to do
    const Utype eps = std::numeric_limits<Utype>::epsilon();

    //size of the polynomial is 0
    if (poly.degree() < 1){
      return;
    }

    bool found = false;                                       // Bool to know if the root was found
    U unpolishRoot;

    U x0 = x;
    U x1 = x0 + std::complex<Utype>(0.1);
    U x2 = x1 + std::complex<Utype>(0.1);

    U f0 = bmt::evaluate_polynomial(&poly.data()[0], x0, poly.size());
    U f1 = bmt::evaluate_polynomial(&poly.data()[0], x1, poly.size());

    if (std::abs(f0) <= eps){                                // Evaluates if the first and second numbers are the root
      unpolishRoot = x0;
      found = true;
    } else if(std::abs(f1) <= eps) {
      unpolishRoot = x1;
      found = true;
    }

    for(int i = 0; i < maxi; ++i){
      if(!(found)){
        U f2 = bmt::evaluate_polynomial(&poly.data()[0], x2, poly.size());

        if(std::abs(f2) <= eps){
          unpolishRoot = x2;
          found = true;
          break;
        }
        // Variables to resolve the muller's method
        U h0 = x1 - x0;
        U h1 = x2 - x1;

        U d0 = (f1 - f0) / h0;
        U d1 = (f2 - f1) / h1;

        U a = (d1 - d0) / (h1 - h0);
        U b = a * h1 + d1;
        U c = f2;

        // Because the method produces two roots, chooses the one with the same sign as b.
        U x3 = (b.real() >= std::complex<Utype>(0).real()) ? x2 - ((U(-2) * c) / (b + sqrt(b * b - U(4) * a * c)))
                                                           : x2 - ((U(-2) * c) / (b - sqrt(b * b - U(4) * a * c)));

        x2 = x3;
        x1 = x2;
        x0 = x1;

        f1 = f2;
        f0 = f1;
      } else {
        break;
      }
    }

    //Polishes if necessary
    if(found && (polish == PolishRoots)){
      anpi::mullerRoots(originalPoly, originalPoly, unpolishRoot, root);
    } else if (found){
      root = unpolishRoot;
    }
  }

  /**
   * Compute the roots of the given polynomial using the Muller method.
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] roots all roots found
   * @param[in] start initial point for finding the roots
   * @param[in] polish indicate if polishing is needed or not.
   *
   * @return the number of roots found
   */
  template<class T,class U>
  void muller(const bmt::polynomial<T>& poly,
              std::vector<U>& roots,
              const PolishEnum polish = DoNotPolish,
              const U start = U()) {

    static_assert(std::is_floating_point<T>::value ||
                  boost::is_complex<T>::value,
                  "T must be floating point or complex");
    static_assert(std::is_floating_point<U>::value ||
                  boost::is_complex<U>::value,
                  "U must be floating point or complex");
    
    bmt::polynomial<T> predeflated = poly;                  // Variable to hold the deflated polynomial
    bmt::polynomial<T> postdeflated;                        // Variable to save the residual when deflated
    
    for(unsigned int i = 0; i < poly.degree(); ++i){

      U residual;
      roots.push_back(start);
      anpi::mullerRoots(predeflated, poly, start, roots.at(i), polish);
      postdeflated = anpi::deflate(predeflated, roots.at(i), residual);   // Deflates the polynomial to find new roots
      predeflated = postdeflated;                                             // The residual is assign as the new polynomial.
    }
  }
}

#endif
