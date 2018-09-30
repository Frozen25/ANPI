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
	void backwardSubs (anpi::Matrix<T>& A,
	                    std::vector <T>& x,
						std::vector <T>& b){
	  size_t i,j;
	  size_t n = A.cols();

	  x = std::vector<T> (n,T(0));
	  T Sum = T(0);
	  
	  x[n-1] = b[n-1]/(A[n-1][n-1]);


	  for (size_t k = (n-1); k > 0 ; --k){
	  	i = k-1;
	  	Sum = b[i];

	  	for(j = i+1 ; j<n ; ++j){
	  		Sum -= A[i][j] * x[j];
	  	}
	  	x[i] = Sum / (A[i][i]);
	  		  	
	  }

	 
	} //backwards subs



	template<typename T>
	void forwardSubs( anpi::Matrix<T>& A,
										std::vector <T>& x,
										std::vector <T>&b){

	  size_t i,j;
	  size_t n = A.cols();

	  T s;

	  for(i = 0; i < n; ++i) {
	        s = T(0);
	        for(j = 0; j < i; ++j) {	                      
	            s = s + A[i][j] * x[j];
	        }
	        x[i] = ( b[i] - s) / A[i][i];
	   }

	}//forwardSubs


	template<typename T>
	void solveLU(const anpi::Matrix<T>& A,
				std::vector <T>& x,
				const std::vector <T>&b){

			Matrix<T> LU;
			size_t size = A.rows();
			LU.allocate(size,size);
			std::vector<size_t> permut;


			lu(A, LU, permut);

			matrix_show(LU);
			std::cout << "permut: ";
			vector_show(permut);

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
			backwardSubs(U, x, Y);

	}


	template<typename T>
	void invert(const Matrix<T>& A,
               Matrix<T>& Ai) {

		const size_t size =  A.cols();
		//anpi::Matrix<T> I = anpi::identityMatrix(size,size);

		std::vector<T> I_j (size,0);
		std::vector<T> AiT_j (size,0);
		anpi::Matrix<T> AiT;

		size_t i,j;
		AiT = A;

		for( i = 0; i< size; ++i){
			I_j[i] = T(1);
			solveLU(A, AiT_j, I_j);
				for( j = 0; j<size; ++j){
					AiT[j][i] = AiT_j[j];
				}
			I_j[i] = T(0);			

		}
	}

	

}//namespace ANPI

#endif