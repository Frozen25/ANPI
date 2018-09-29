
#include <boost/test/unit_test.hpp>

#include "LUCrout.hpp"
#include "LUDoolittle.hpp"
#include "Solver.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
    namespace test {


        template<typename T>
        void invertTest(const std::function<void(const Matrix<T>& A,
                                                 Matrix<T>& Ai)>& invert) {    //Test the unpack doolittle method

            //factorized matrix in LU form with Doolittle
            Matrix<T> AA;

            Matrix<T> AAi;



            {
                //Test with 3x3 matrix
                AA = { {1, 2, 3},{0,1,4},{5,6,0} };
                invert(AA, AAi);
                //Real Ai matrix
                Matrix<T> Ai_real = {{-24, 18, 5},{20, -15, -4},{-5, 4, 1}};
                
                //Test each element one by one
                for (size_t i=0;i<AA.rows();++i) {
                    for (size_t j=0;j<AA.cols();++j) {
                        BOOST_CHECK(AAi(i,j)==Ai_real(i,j));
                    }
                }
            }
        }


    } // test
}  // anpi






BOOST_AUTO_TEST_SUITE( inverter )

    BOOST_AUTO_TEST_CASE(inverters)
    {
        anpi::test::invertTest<float>(anpi::invert<float>);

        anpi::test::invertTest<double>(anpi::invert<double>);

    }


BOOST_AUTO_TEST_SUITE_END()
