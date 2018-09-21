/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_JENKINS_TRAUB_HPP
#define ANPI_JENKINS_TRAUB_HPP

#include <vector>
#include <type_traits>
#include <cctype>
#include <cmath>
#include <cfloat>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>


#define MAXDEGREE 100
#define MDP1 MAXDEGREE+1

void rpoly_ak1(T op[MDP1], int* Degree, T zeror[MAXDEGREE], T zeroi[MAXDEGREE]);
void Fxshfr_ak1(int L2, int* NZ, T sr, T bnd, T K[MDP1], int N, T p[MDP1], int NN, T qp[MDP1], T* lzi, T* lzr, T* szi, T* szr);
void QuadSD_ak1(int NN, T u, T v, T p[MDP1], T q[MDP1], T* a, T* b);
int calcSC_ak1(int N, T a, T b, T* a1, T* a3, T* a7, T* c, T* d, T* e, T* f, T* g, T* h, T K[MDP1], T u, T v, T qk[MDP1]);
void nextK_ak1(int N, int tFlag, T a, T b, T a1, T* a3, T* a7, T K[MDP1], T qk[MDP1], T qp[MDP1]);
void newest_ak1(int tFlag, T* uu, T* vv, T a, T a1, T a3, T a7, T b, T c, T d, T f, T g, T h, T u, T v, T K[MDP1], int N, T p[MDP1]);
void QuadIT_ak1(int N, int* NZ, T uu, T vv, T* szr, T* szi, T* lzr, T* lzi, T qp[MDP1], int NN, T* a, T* b, T p[MDP1], T qk[MDP1], T* a1, T* a3, T* a7, T* d, T* e, T* f, T* g, T* h, T K[MDP1]);
void RealIT_ak1(int* iFlag, int* NZ, T* sss, int N, T p[MDP1], int NN, T qp[MDP1], T* szr, T* szi, T K[MDP1], T qk[MDP1]);
void Quad_ak1(T a, T b1, T c, T* sr, T* si, T* lr, T* li);


namespace anpi {

  namespace bmt=boost::math::tools; // for polynomial
  
  /**
   * Compute the roots of the given polynomial using the Jenkins-Traub method.
   * @param[in] poly polynomial to be analyzed for roots
   * @param[out] roots all roots found
   * @param[in] start initial point for finding the roots
   * @param[in] polish indicate if polishing is needed or not.
   *
   * @return the number of roots found
   */
  template <class T, class U>
void jenkinsTraub(const bmt::polynomial<T>& poly,
                    std::vector<U>& roots) {
    
    static_assert(std::is_floating_point<T>::value ||
                  boost::is_complex<T>::value,
                  "T must be floating point or complex");
    static_assert(std::is_floating_point<U>::value ||
                  boost::is_complex<U>::value,
                  "U must be floating point or complex");


    /*
    Código no funcional para reales: JenkinsR.hpp
    Código no funcional para reales y complejos: JenkinsRC.hpp



    
   int degree = poly.degree();
   T polynom[degree];
   for (int i = 0; i <= degree ;++i ){
     polynom[i] = poly[i];
   }
   ZeroReals[degree];
   ZeroImag [degree];

   rpoly_ak1( polynom, &degree, ZeroReals, ZeroImag );

   std::vector<T> vectR(ZeroReals, ZeroReals + sizeof ZeroReals / sizeof ZeroReals[0]);
   std::vector<T> vectI(ZeroImag,  ZeroImag +  sizeof ZeroImag /  sizeof ZeroImag[0]);



    roots.insert( roots.end(), vectR.begin(), vectR.end() );
    roots.insert( roots.end(), vectI.begin(), vectI.end() );

    */

    throw Exception("Not implemented yet!");
  }



}//end namespace anpi



#endif
