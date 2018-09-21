#include <cmath>
#include <limits>
#include <functional>
#include <vector>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "LU.hpp"
#include "Utilities.hpp"

#ifndef ANPI_SOLVER_HPP
#define ANPI_SOLVER_HPP

namespace anpi {

	template<typename T>
	void backwardsSubs (anpi::Matrix<T>& A,
	                    std::vector <T>& x,
						std::vector <T>& b){
	  T i,j;
	  size_t n = A.cols();

	  for (i = 0; i<n ; ++i){
		x[i] = b[n-1]/A[n-1][n-1];
	  }

	  T Sum;

	  for (i=n-1; i>=0; --i){
	    Sum = T(0);

	    for (j=i+1; j<n; ++j){
	      Sum -= A[i][j]*x[j];
	    }    

	      x[i] = Sum / A[i][i];
    }
	} //backwards subs



	template<typename T>
	void forwardSubs( anpi::Matrix<T>& A,
										std::vector <T>& x,
										std::vector <T>&b){

	  size_t i,j;
	  size_t n = A.cols();

	  T s;

	  for(i = 0; i < n; i++) {
	        s = T(0);
	        for(j = 0; j < i; j++) {	                      
	            s = s + a[i][j] * x[j];
	        }
	        x[i] = ( b[i] - s) / a[i][i];
	   }

	}


	template<typename T>
	void solveLU(const anpi::Matrix<T>& A,
				std::vector <T>& x,
				const std::vector <T>&b){

			Matrix<T> LU;
			size_t size = A.rows();
			LU.allocate(size,size);
			std::vector<size_t> permut;

			lu(A, LU, permut);

			Matrix<T> L, U;
			L.allocate(size, size);
			U.allocate(size, size);
			unpack(LU, L, U);



			std::vector<T> Y (size,0);

			std::vector<T> Bperm (size,0);
			for(size_t i = 0; i<size; ++i){
				Bperm[i] = b[permut[i]];
			}

			forwardSubs(L, Y, Bperm);
			backwardsSubs(U, x, Y);

	}


	template<typename T>
	void invert(const Matrix<T>& A,
               Matrix<T>& Ai) {

		const size_t size = (size_t) A.cols();
		anpi::Matrix<T> I = anpi::identityMatrix(size,size);

		std::vector<T> I_j (size,0);
		std::vector<T> AiT_j (size,0);
		anpi::Matrix<T> AiT;

		AiT = A;

		for(size_t i = 0; i< size; ++i){
			I_j[i] = T(1);
			solveLU(A, AiT_j, I_j);
				for(size_t j = 0; j<size; ++j){
					Matrix[j][i] = AiT_j[j];
				}
			I_j[i] = T(0);			

		}
	}

	/**
	 * Function that obtains the values of x with a given vector b and a matrix A
	 * @tparam T type of data
	 * @param[in] A base matrix that contains the equations
	 * @param[out] x result of the equations
	 * @param[in] b factors of the matrix
	 * @return
	 */
    template<typename T>
    bool solveQR (const anpi::Matrix<T>& A,
                  std::vector<T>& x,
                  const std::vector<T>& b) {

        Matrix<T> Q;
        Matrix<T> R;

        anpi::qr(A, Q, R);

        Matrix<T> QT = Q;
        QT.transpose();

        std::vector<T> b_primed = QT*b;
        size_t last_position = b.size()-1;

        x.resize(last_position+1);
        x[last_position] = b_primed[last_position]/R(last_position,last_position);

        for (size_t i = last_position-1; i >= 0; --i) {
            T sum;
            for (size_t j = last_position; j > i; --j) {
                sum += R(i,j)*x[j];
            }
            x[i] = (b_primed[i]-sum)/R(i,i);

            if (i == 0){
                return true;
            }
        }

        return false;
    }

}//namespace ANPI

#endif