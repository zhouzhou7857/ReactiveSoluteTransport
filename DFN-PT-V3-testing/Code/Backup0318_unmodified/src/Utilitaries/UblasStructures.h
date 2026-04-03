/*
 * UblasStructures.h
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */

#ifndef UBLASSTRUCTURES_H_
#define UBLASSTRUCTURES_H_

#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::column_major, 0,boost::numeric::ublas::unbounded_array<int>, boost::numeric::ublas::unbounded_array<double> > ublas_matrix;
typedef boost::numeric::ublas::vector<double> ublas_vector;

void print_matrix(ublas_matrix,int,int,int option=0);
void print_vector(ublas_vector,int);


#endif /* UBLASSTRUCTURES_H_ */
