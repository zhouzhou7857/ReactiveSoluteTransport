/*
 * LinearSystem.cpp
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */

#include "LinearSystem.h"

#include <iostream>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace ublas = boost::numeric::ublas;
namespace umf = boost::numeric::bindings::umfpack;
using namespace std;

// return the solution x of the linear system Ax = b (A = mat_syst, b = vect_syst and x = vect_sol)
ublas_vector LinearSystemSolving(ublas_matrix mat_syst, ublas_vector vect_syst){

	// 1. Variables definition
	int Ny = vect_syst.size();
	ublas::vector<double> X(Ny);

	// function implemented for 2D matrix
	if (mat_syst.size2()==vect_syst.size()){
		// 2. System solving
		umf::symbolic_type<double> Symbolic;
		umf::numeric_type<double> Numeric;
		umf::symbolic (mat_syst, Symbolic);
		umf::numeric (mat_syst, Symbolic, Numeric);
		umf::solve (mat_syst, X, vect_syst, Numeric);

	}
	else{cout << "WARNING in LinearSystemSolving (LinearSystemUtilitaries.cpp): system not well defined" << endl;}

	return X;
}






