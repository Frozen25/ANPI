/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Andres Ramirez-Quiros
 * @Date  : 17.09.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>

#include "Utilities.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_QR_HPP
#define ANPI_QR_HPP

namespace anpi {

   /**
   * Function no get the norm of a vector
   * @tparam T type of data
   * @param a vector to get the norm
   * @return the norm fo the vector
   */
  template<typename T>
  T norm(const std::vector<T>& a) {
    T sum = T(0);
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * a[i];			// Access a vector in i
    }
    return sqrt(sum);
  }

  /**
   * Function that scalates a vector with a given value
   * @tparam T type of data
   * @param[in] factor number by which scalate the function
   * @param[in,out] a vector to scalate
   */
  template<typename T>
  void rescale(T factor, std::vector<T>& a) {
    for (size_t i = 0; i < a.size(); ++i) a[i] /= factor;
  }

  /**
   * Function that scales its vector according to its norm
   * @tparam T type of data
   * @param[in,out] a vector scaled
   */
  template<typename T>
  void rescale_unit(std::vector<T>& a) {
    T factor = norm(a);
    rescale(factor, a);
  }

  /**
   * Function that does c = a + b * s, necessary for the QR calculation
   * @tparam T type of data
   * @param[in] a vector a
   * @param[in] b vector b
   * @param[in] s vector s
   * @param[out] c vector c where the result is stored
   */
  template<typename T>
  void vmadd(const std::vector<T>& a,
             const std::vector<T>& b, T s, std::vector<T>& c) {

    if (c.size() != a.size() or c.size() != b.size()) {
        throw anpi::Exception("[vmadd]: vector sizes don't match\n");
    }

    for (size_t i = 0; i < c.size(); ++i)
        c[i] = a[i] + s * b[i];
  }

  /**
   * Function that computes de HouseHolder Factor using I - 2*v*v^T
   * @tparam T type of data
   * @tparam Alloc
   * @param[out] a Matrix to store the householder matrix
   * @param[in] v vector necessary for the formula
   */
  template<typename T,class Alloc>
  void compute_householder_factor(Matrix<T,Alloc>& a,
                                  const std::vector<T>& v) {

    size_t n = (size_t) v.size();
    a.allocate(n,n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j)
        a(i,j) = -2 *  v[i] * v[j];
    }

    for (size_t i = 0; i < n; i++)
      a(i,i) += 1;
  }

  /**
   * Function that decomposes a matrix A to a Q matrix (orthogonal) and a R matrix (triangular superior)
   * @tparam T type of data
   * @param[in] A base matrix A
   * @param[out] Q matrix Q
   * @param[out] R matrix R
   */
  template<typename T>
  void qr(const anpi::Matrix<T>& A,
          anpi::Matrix<T>& Q,
          anpi::Matrix<T>& R ) {

    if (A.rows() != A.cols()) throw anpi::Exception("Matrix is not a square!");

    size_t m = A.rows();
    size_t n = A.cols();

    // array of factor Q1, Q2, ... Qm
    std::vector<anpi::Matrix<T>> qv(m);

    // temp array
    anpi::Matrix<T> z(A);
    anpi::Matrix<T> z1;

    for (size_t k = 0; k < n && k < m-1; k++) {

      std::vector<T> e(m), x(m);
      T a;

      //compute minor
      z1.compute_minor(z, k);

      // extract k-th column into x
      //z1.extract_column(x, k);
      x = z1.column(k);

      a = norm(x);
      if (A(k,k) > 0) {
        a = -a;
      }

      for(size_t i = 0; i < e.size(); ++i) {
        e[i] = (i == k) ? 1 : 0;
      }

      // e = x + a*e
      vmadd(x, e, a, e);

      // e = e / ||e||
      rescale_unit(e);

      // qv[k] = I - 2 *e*e^T
      compute_householder_factor(qv[k], e);

      // z = qv[k] * z1
      z = qv[k] * z1;
    }

    Q = qv[0];

    // after this loop, we will obtain Q (up to a transpose operation)
 	for (size_t i = 1; i < n && i < m - 1; i++) {
 
      z1 = qv[i] * Q;
      Q = z1;
 
  	}
 
  	R = Q * A;
  	Q.transpose();

  }
}

#endif //ANPI_QR_HPP
