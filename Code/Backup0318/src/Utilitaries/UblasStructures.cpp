/*
 * UblasStructures.cpp
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */

#include "UblasStructures.h"

using namespace std;

void print_matrix(ublas_matrix ublas_mat, int Nx, int Ny,int option){
	cout << "Matrix values = " << endl;
	if (option==0){
		for (int i=0;i<Nx;i++){
			for (int j=0;j<Ny;j++){
				cout << ublas_mat(i,j) << " ";
			}
			cout << endl;
		}
	}
	else if (option==1){
		for (int j=Ny-1;j>=0;j--){
			for (int i=0;i<Nx;i++){
				cout << ublas_mat(i,j) << " ";
			}
			cout << endl;
		}
	}
	else{cout << "WARNING in print_matrix (UblasStructures.cpp): option not defined" << endl;}
}

void print_vector(ublas_vector ublas_vect, int N){
	cout << "Vector values = " << endl;
	for (int i=0;i<N;i++){
		cout << ublas_vect(i) << endl;
	}
}
