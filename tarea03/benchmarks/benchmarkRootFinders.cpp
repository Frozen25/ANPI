/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @post editor: Andres Ramirez Quiros
 * @date   27.08.2018
 */
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <PlotPy.hpp>

#include "Exception.hpp"

/**
 * Load the root finders themselves
 */
#include "RootBisection.hpp"
#include "RootInterpolation.hpp"
#include "RootSecant.hpp"
#include "RootBrent.hpp"
#include "RootNewtonRaphson.hpp"
#include "RootRidder.hpp"

#include "Allocator.hpp"

namespace anpi {
  // namespace encapsulating benchmarking functionality
  namespace bm { 

    /// Square of a number
    template<typename T>
    inline T sqr(const T x) { return x*x; }

    /// Cube of a number
    template<typename T>
    inline T cube(const T x) { return x*x*x; }
    
    /// First testing function for roots |x|=e^(-x)
    template<typename T>
    T t1(const T x)  { return std::abs(x)-std::exp(-x); }

    /// Second testing function for roots e^(-x²) = e^(-(x-3)²/3 )
    template<typename T>
    T t2(const T x) { return std::exp(-x*x) - std::exp(-sqr(x-T(3))/T(3)); }

    /// Third testing function for roots x² = atan(x)
    template<typename T>
    T t3(const T x)  { return x*x-std::atan(x); }

      /// Fourth testing function for roots (x-2)^3 +0.01(x-2)
    template<typename T>
    T t4(const T x)  { const T x0=x-T(2); return cube(x0) + T(0.01)*x0; }

    /**
     * Wrapper class to count function calls
     *
     * This wrapper fulfills the requirements to act as a
     * std::function<T(T)>, and it simply counts the number
     * of calls made to the operator(), before calling
     * the functor provided at construction time.
     */
    template<typename T>
    class CallCounter {
    protected:
      /// Maximum allowed size for the square matrices
      mutable size_t _counter;
      
      std::function<T(T)> _f;
    public:
      /// Construct
      CallCounter(std::function<T(T)> f) : _counter(0u),_f(f) {}
      
      /// Access the counter
      inline size_t counter() const {return _counter;}
      
      /// Reset the counter
      inline void reset() { _counter = 0u; }
      
      /// Call the function
      T operator()(const T x) {
        ++_counter;
        return _f(x);
      }
      
      /// Call the function
      T operator()(const T x) const { ++_counter; return _f(x); }
    };

    /**
     * Test the given _closed_ root finder
     *
     * The solver must be itself a std::function expecting another
     * std::function (the one whose roots are being looked for), the
     * two limits of the interval enclosing the root, and the
     * tolerance.
     *
     * The tolerances will start from "start", then progressing with
     *   eps = eps*factor
     * until the end value is reached.
     */
    template<typename T>
    void rootBench(const std::function<T(const std::function<T(T)>&,
                                         T,
                                         T,
                                         const T)>& solver,
                   const T start,
                   const T end,
                   const T factor) {

      if ( (factor >= static_cast<T>(1)) &&
           (factor < static_cast<T>(0)) ) {
        throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
      }

      // Alias of the function type, for which the roots are being looked for.
      typedef std::function<T(T)> f_type;

      // Try a series of tolerances
      for (T eps=start; eps>end; eps*=factor) {
        std::cout << "eps=" << eps << "; ";

        // Create an std::function instance, which wraps the function
        // t1 with the function counter
        f_type c1(CallCounter<T>(t1<T>));
        solver(c1,T(0),T(2),eps);
        std::cout << c1.template target< CallCounter<T> >()->counter() << "; ";

        // now the same with function t2
        f_type c2(CallCounter<T>(t2<T>));
        solver(c2,T(0),T(2),eps);
        std::cout << c2.template target< CallCounter<T> >()->counter() << "; ";

        // now the same with function t3
        f_type c3(CallCounter<T>(t3<T>));
        solver(c3,T(0),T(0.5),eps);
        std::cout << c3.template target< CallCounter<T> >()->counter() << "; ";

        // now the same with function t4
        f_type c4(CallCounter<T>(t4<T>));
        solver(c4,T(1),T(3),eps);
        std::cout << c4.template target< CallCounter<T> >()->counter() << std::endl;
      }
    }

    /**
     * Test the given _open_ root finder
     *
     * The solver must be itself a std::function expecting another
     * std::function (the one whose roots are being looked for), the
     * starting root guess, and the tolerance.
     */
    template<typename T>
    void rootBench(const std::function<T(const std::function<T(T)>&,
                                         T,
                                         const T)>& solver,
                   const T start,
                   const T end,
                   const T factor) {

      if ( (factor >= static_cast<T>(1)) &&
           (factor < static_cast<T>(0)) ) {
        throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
      }
     
      // Alias of the function type, for which the roots are being looked for.
      typedef std::function<T(T)> f_type;
      
      // Try a series of tolerances
      for (T eps=start; eps>end; eps*=factor) {
        std::cout << "eps=" << eps << "; ";

        // Create an std::function instance, which wraps the function
        // t1 with the function counter       
        f_type c1(CallCounter<T>(t1<T>));
        solver(c1,T(0),eps);
        std::cout << c1.template target< CallCounter<T> >()->counter() << "; ";

        // now the same with function t2
        f_type c2(CallCounter<T>(t2<T>));
        solver(c2,T(2),eps);
        std::cout << c2.template target< CallCounter<T> >()->counter() << "; ";
        
        // now the same with function t3
        f_type c3(CallCounter<T>(t3<T>));
        solver(c3,T(0),eps);
        std::cout << c3.template target< CallCounter<T> >()->counter() << "; ";
        
        // now the same with function t4
        f_type c4(CallCounter<T>(t4<T>));
        solver(c4,T(1),eps);
        std::cout << c4.template target< CallCounter<T> >()->counter() << std::endl;
      }
    }

    /**
     * Method that gets the number of calls a method did by trying to find the zero of a function
     * @tparam T: type of the parameter, normally float or double
     * @param solver: method to execute with the tests
     * @param start: starting point of the error interval
     * @param end: ending point of the error interval
     * @param factor: step by which the number is reduce in each iteration
     * @param func: number of the function to execute. Ranges from 0 to 3. See the available functions for this.
     * @param ans: pointer to the vector where the answer is going to be stored
     */
    template<typename T>
    void oneRoot(const std::function<T(const std::function<T(T)>&,
                                         T,
                                         const T)>& solver,
                   const T start,
                   const T end,
                   const T factor,
                   const int func,
                   std::vector<double>& ans) {

        if ( (factor >= static_cast<T>(1)) &&
             (factor < static_cast<T>(0)) ) {
            throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
        }

        // Alias of the function type, for which the roots are being looked for.
        typedef std::function<T(T)> f_type;

        // counter to access the vector
        int i = 0;

        // Try a series of tolerances
        for (T eps=start; eps>end; eps*=factor) {
            std::cout << "eps=" << eps << "; ";

            switch (func) {
                case 0: {
                    // Create an std::function instance, which wraps the function
                    // t1 with the function counter
                    f_type c(CallCounter<T>(t1<T>));
                    solver(c,T(0),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                case 1: {
                    // now the same with function t2
                    f_type c(CallCounter<T>(t2<T>));
                    solver(c,T(2),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                case 2: {
                    // now the same with function t3
                    f_type c(CallCounter<T>(t3<T>));
                    solver(c,T(0),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                case 3: {
                    // now the same with function t4
                    f_type c(CallCounter<T>(t4<T>));
                    solver(c,T(1),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                default:
                    break;
            }

            ++i;
        }
    }

   /**
    * Method that gets the number of calls a method did by trying to find the zero of a function
    * @tparam T: type of the parameter, normally float or double
    * @param solver: method to execute with the tests
    * @param start: starting point of the error interval
    * @param end: ending point of the error interval
    * @param factor: step by which the number is reduce in each iteration
    * @param func: number of the function to execute. Ranges from 0 to 3. See the available functions for this.
    * @param ans: pointer to the vector where the answer is going to be stored
    */
    template<typename T>
    void oneRoot(const std::function<T(const std::function<T(T)>&,
                                                 T,
                                                 T,
                                                 const T)>& solver,
                   const T start,
                   const T end,
                   const T factor,
                   const int func,
                   std::vector<double>& ans) {

        if ( (factor >= static_cast<T>(1)) &&
             (factor < static_cast<T>(0)) ) {
            throw anpi::Exception("Invalid factor.  It must be between 0 and 1");
        }

        // Alias of the function type, for which the roots are being looked for.
        typedef std::function<T(T)> f_type;

        // counter to access the vector
        int i = 0;

        // Try a series of tolerances
        for (T eps=start; eps>end; eps*=factor) {
            std::cout << "eps=" << eps << "; ";

            switch (func) {
                case 0: {
                    // Create an std::function instance, which wraps the function
                    // t1 with the function counter
                    f_type c(CallCounter<T>(t1<T>));
                    solver(c,T(0),T(2),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                case 1: {
                    // now the same with function t2
                    f_type c(CallCounter<T>(t2<T>));
                    solver(c,T(0),T(2),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                case 2: {
                    // now the same with function t3
                    f_type c(CallCounter<T>(t3<T>));
                    solver(c,T(0),T(0.5),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                case 3: {
                    // now the same with function t4
                    f_type c(CallCounter<T>(t4<T>));
                    solver(c,T(1),T(3),eps);
                    std::cout << c.template target< CallCounter<T> >()->counter() << std::endl;
                    ans[i] = c.template target< CallCounter<T> >()->counter();
                    break;
                }
                default:
                    break;
            }

            ++i;
        }
    }

    /**
     * Benchmark all solvers using a range of tolerances geometrically changing
     * multiplying from the start point until the end with the given factor
     */
    template<typename T>
    void allSolvers(const T start,const T end,const T factor) {

      std::cout << "Bisection" << std::endl;
      anpi::bm::rootBench<T>(anpi::rootBisection<T>,start,end,factor);

      std::cout << "Interpolation" << std::endl;
      anpi::bm::rootBench<T>(anpi::rootInterpolation<T>,start,end,factor);

      std::cout << "Secant" << std::endl;
      anpi::bm::rootBench<T>(anpi::rootSecant<T>,start,end,factor);

      std::cout << "NewtonRaphson" << std::endl;
      anpi::bm::rootBench<T>(anpi::rootNewtonRaphson<T>,start,end,factor);

      std::cout << "Brent" << std::endl;
      anpi::bm::rootBench<T>(anpi::rootBrent<T>,start,end,factor);

      std::cout << "Ridder" << std::endl;
      anpi::bm::rootBench<T>(anpi::rootRidder<T>,start,end,factor);
    }

    /**
     * Method that gets the number of calls a method did by trying to find the zero of a function
     * @tparam T: type of the parameter, normally float or double
     * @param start: starting point of the error interval
     * @param end: ending point of the error interval
     * @param factor: step by which the number is reduce in each iteration
     * @param method: method to use to resolve a functions. Ranges from 0 to 5. See the available methods.
     * @param func: number of the function to execute. Ranges from 0 to 3. See the available functions for this.
     * @param ans: pointer to the vector where the answer is going to be stored
     */
    template<typename T>
    void oneSolver(const T start,const T end,const T factor,
            const int method, const int func, std::vector<double>& ans) {

      switch (method){
          case 0: {
              std::cout << "Bisection" << std::endl;
              anpi::bm::oneRoot<T>(anpi::rootBisection<T>,start,end,factor,func,ans);
              break;
          }
          case 1: {
              std::cout << "Interpolation" << std::endl;
              anpi::bm::oneRoot<T>(anpi::rootInterpolation<T>,start,end,factor,func,ans);
              break;
          }
          case 2: {
              std::cout << "Secant" << std::endl;
              anpi::bm::oneRoot<T>(anpi::rootSecant<T>,start,end,factor,func,ans);
              break;
          }
          case 3: {
              std::cout << "NewtonRaphson" << std::endl;
              anpi::bm::oneRoot<T>(anpi::rootNewtonRaphson<T>,start,end,factor,func,ans);
              break;
          }
          case 4: {
              std::cout << "Brent" << std::endl;
              anpi::bm::oneRoot<T>(anpi::rootBrent<T>,start,end,factor,func,ans);
              break;
          }
          case 5: {
              std::cout << "Ridder" << std::endl;
              anpi::bm::oneRoot<T>(anpi::rootRidder<T>,start,end,factor,func,ans);
              break;
          }
          default:
              break;
      }
    }
  } // bm
}  // anpi

BOOST_AUTO_TEST_SUITE( RootFinders )

/**
 * Instantiate and test the methods of the Matrix class
 */
/* Test not necessary because the ones under execute the same when buildin the graph
BOOST_AUTO_TEST_CASE( RootFinders ) {

  // Benchmark the solvers using float
  std::cout << "<float>" << std::endl;
  anpi::bm::allSolvers<float>(0.1f,1.e-7f,0.125f);
  
  // Benchmark the solvers using double
  std::cout << "<double>" << std::endl;
  anpi::bm::allSolvers<double>(0.1f,1.e-15f,0.125f);
}
*/
/////////////////////////////////////////////////////////////// Graph of each function with all the methods in float
BOOST_AUTO_TEST_CASE( Function_1_float ) {
    static anpi::Plot2d<double> plotter;
    plotter.initialize(1);

    const float start = 0.1f;
    const float end = 1.e-7f;
    const float factor = 0.125f;
    const int maxi = 7;

    std::vector<double> x_axis(maxi), y_axis(maxi);
    x_axis[0] = start;

    for (int i = 1; i < maxi; ++i) {
        x_axis[i] = x_axis[i-1] * factor;
    }

    ::anpi::bm::oneSolver(start,end,factor,0,0,y_axis);
    plotter.plot(x_axis,y_axis,"Biseccion","red");

    ::anpi::bm::oneSolver(start,end,factor,1,0,y_axis);
    plotter.plot(x_axis,y_axis,"Interpolacion","blue");

    ::anpi::bm::oneSolver(start,end,factor,2,0,y_axis);
    plotter.plot(x_axis,y_axis,"Secante","green");

    ::anpi::bm::oneSolver(start,end,factor,3,0,y_axis);
    plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

    ::anpi::bm::oneSolver(start,end,factor,4,0,y_axis);
    plotter.plot(x_axis,y_axis,"Brent","yellow");

    ::anpi::bm::oneSolver(start,end,factor,5,0,y_axis);
    plotter.plot(x_axis,y_axis,"Ridder","purple");

    plotter.setXLabel("Tolerancia de respuesta");
    plotter.setYLabel("Numero de iteraciones");
    plotter.setTitle("Solucion de la funcion 1 utilizando diferentes metodos con float");

    plotter.show();
}

BOOST_AUTO_TEST_CASE( Function_2_float ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const float start = 0.1f;
        const float end = 1.e-7f;
        const float factor = 0.125f;
        const int maxi = 7;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,1,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,1,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,1,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,1,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,1,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,1,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 2 utilizando diferentes metodos con float");

        plotter.show();
}

BOOST_AUTO_TEST_CASE( Function_3_float ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const float start = 0.1f;
        const float end = 1.e-7f;
        const float factor = 0.125f;
        const int maxi = 7;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,2,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,2,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,2,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,2,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,2,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,2,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 3 utilizando diferentes metodos con float");

        plotter.show();
}

BOOST_AUTO_TEST_CASE( Function_4_float ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const float start = 0.1f;
        const float end = 1.e-7f;
        const float factor = 0.125f;
        const int maxi = 7;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,3,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,3,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,3,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,3,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,3,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,3,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 4 utilizando diferentes metodos con float");

        plotter.show();
}

/////////////////////////////////////////////////////////////// Graph of each function with all the methods in double
BOOST_AUTO_TEST_CASE( Function_1_double ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const double start = 0.1f;
        const double end = 1.e-15f;
        const double factor = 0.125f;
        const int maxi = 16;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,0,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,0,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,0,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,0,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,0,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,0,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 1 utilizando diferentes metodos con double");

        plotter.show();
}

BOOST_AUTO_TEST_CASE( Function_2_double ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const double start = 0.1f;
        const double end = 1.e-15f;
        const double factor = 0.125f;
        const int maxi = 16;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,1,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,1,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,1,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,1,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,1,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,1,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 2 utilizando diferentes metodos con double");

        plotter.show();
}

BOOST_AUTO_TEST_CASE( Function_3_double ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const double start = 0.1f;
        const double end = 1.e-15f;
        const double factor = 0.125f;
        const int maxi = 16;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,2,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,2,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,2,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,2,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,2,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,2,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 3 utilizando diferentes metodos con double");

        plotter.show();
}

BOOST_AUTO_TEST_CASE( Function_4_double ) {
        static anpi::Plot2d<double> plotter;
        plotter.initialize(1);

        const double start = 0.1f;
        const double end = 1.e-15f;
        const double factor = 0.125f;
        const int maxi = 16;

        std::vector<double> x_axis(maxi), y_axis(maxi);
        x_axis[0] = start;

        for (int i = 1; i < maxi; ++i) {
            x_axis[i] = x_axis[i-1] * factor;
        }

        ::anpi::bm::oneSolver(start,end,factor,0,3,y_axis);
        plotter.plot(x_axis,y_axis,"Biseccion","red");

        ::anpi::bm::oneSolver(start,end,factor,1,3,y_axis);
        plotter.plot(x_axis,y_axis,"Interpolacion","blue");

        ::anpi::bm::oneSolver(start,end,factor,2,3,y_axis);
        plotter.plot(x_axis,y_axis,"Secante","green");

        ::anpi::bm::oneSolver(start,end,factor,3,3,y_axis);
        plotter.plot(x_axis,y_axis,"Newton-Raphson","black");

        ::anpi::bm::oneSolver(start,end,factor,4,3,y_axis);
        plotter.plot(x_axis,y_axis,"Brent","yellow");

        ::anpi::bm::oneSolver(start,end,factor,5,3,y_axis);
        plotter.plot(x_axis,y_axis,"Ridder","purple");

        plotter.setXLabel("Tolerancia de respuesta");
        plotter.setYLabel("Numero de iteraciones");
        plotter.setTitle("Solucion de la funcion 4 utilizando diferentes metodos con double");

        plotter.show();
}

BOOST_AUTO_TEST_SUITE_END()
