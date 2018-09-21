//
// Created by ubuntu on 18/09/18.
//

/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Andres Ramirez-Quiros
 * @Date  : 18.09.2018
 */

#include <boost/test/unit_test.hpp>

#include "QR.hpp"
#include "Utilities.hpp"
#include "Solver.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
    namespace test {

        /// Test the given closed root finder
        template<typename T>
        void qrTest(const std::function<void(const Matrix<T>&, Matrix<T>&, Matrix<T>&)>& QR,
                                             const Matrix<T>& A) {

            // The result
            Matrix<T> Q;
            Matrix<T> R;
            const T eps = std::numeric_limits<T>::epsilon()*10;

            // Test decomposition
            {
                anpi::qr(A, Q, R);

                Matrix<T> Ar=Q*R;

                BOOST_CHECK(Ar.rows()==A.rows());
                BOOST_CHECK(Ar.cols()==A.cols());

                // Check if A = Ar
                //BOOST_CHECK(A==Ar);
                for (size_t i=0;i<Ar.rows();++i) {
                    for (size_t j=0;j<Ar.cols();++j) {
                        BOOST_CHECK(std::abs(Ar(i,j)-A(i,j)) < eps);
                    }
                }

                Matrix<T> QT=Q;
                QT.transpose();

                Matrix<T> QQT=Q*QT;

                // Check if Q*QT = I
                //BOOST_CHECK(Q*Q.transpose() == anpi::identityMatrix);
                for (size_t i=0;i<QQT.rows();++i) {
                    for (size_t j=0;j<QQT.cols();++j) {
                        if (i==j) {
                          BOOST_CHECK(std::abs(QQT(i,j)-1) < eps);
                        } else {
                          BOOST_CHECK(std::abs(QQT(i,j)) < eps);
                        }
                    }
                }

                // R es triangular superior
                for (size_t i=0;i<R.rows();++i) {
                    for (size_t j=0;j<R.cols();++j) {
                        if (i>j) {
                            BOOST_CHECK(std::abs(R(i,j)) < eps);
                        }
                    }
                }
            }
        }

        template<typename T>
        void qrTestSolver(const std::function<void(const Matrix<T>&, Matrix<T>&, Matrix<T>&)>& QR,
                          const Matrix<T>& A, const std::vector<T>& b) {

            std::vector<T> x;
            anpi::solveQR(A,x,b);
        }

    } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( QRTEST )

    BOOST_AUTO_TEST_CASE(QR)
    {
        const anpi::Matrix<float> Af = { {2,0,1,2},{-1,-2,1,2},{1,1,1,1},{-1,-1,0,1} };
        const anpi::Matrix<double> Ad = { {2,0,1,2},{-1,-2,1,2},{1,1,1,1},{-1,-1,0,1} };

        anpi::test::qrTest<float>(anpi::qr<float>, Af);
        anpi::test::qrTest<double>(anpi::qr<double>, Ad);
    }

    BOOST_AUTO_TEST_CASE(SOLVER)
    {
        const anpi::Matrix<float> Af = { {2,0,1,2},{-1,-2,1,2},{1,1,1,1},{-1,-1,0,1} };
        const std::vector<float> bf = {4,3,3,4};

        const anpi::Matrix<double> Ad = { {2,0,1,2},{-1,-2,1,2},{1,1,1,1},{-1,-1,0,1} };
        const std::vector<double> bd = {4,3,3,4};

        anpi::test::qrTestSolver<float>(anpi::qr<float>, Af, bf);
        anpi::test::qrTestSolver<double>(anpi::qr<double>, Ad, bd);
    }


BOOST_AUTO_TEST_SUITE_END()
