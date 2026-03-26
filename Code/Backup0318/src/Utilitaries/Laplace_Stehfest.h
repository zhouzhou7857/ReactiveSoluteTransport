/*
 * Laplace_Stehfest.h
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#ifndef LAPLACE_STEHFEST_H_
#define LAPLACE_STEHFEST_H_


#include <vector>

/** @addtogroup UTILITARY_MATH UTILITARY MATH Laplace Transform
*  @{
*/

// Stehfest multiplicator coefficients
std::vector<double> Stehfest_Calcul_Vi( int N);

// Sampling computation
std::vector<double> Stehfest_Sampling( int N, double T);

// Inverse Laplace transform by Stehfest algorithm
double Stehfet_Inverse_Laplace( double T, std::vector<double> & fp, int N, std::vector<double> & Vi);

/** Compute Laplace inverse for "formula" at time "t" with success written in "test" */
/** Main function perfroming the time integration*/
template<class formula>
double Laplace_Inverse_Stehfest(formula & f, double t, bool & test){

	// Number of Laplace sampling point for a single time
	int Np = 8;

	// Stehfest coefficients
	std::vector<double> Vi = Stehfest_Calcul_Vi( Np);
	// Results in Laplace space (Hp)
	std::vector<double> Hp(Np+1);
	// Sampling in Laplace space function of p
	std::vector<double> p = Stehfest_Sampling( Np, t);
	// Solving of the Np systems in Laplace space
	for( int j = 1;j <=Np; j++){
		Hp[j] = f.compute_Laplace_p(p[j], test);
		if(!test) return 0;
	}
	// Inversion back in temporal domain
	double h = Stehfet_Inverse_Laplace(t, Hp, Np, Vi);
	return h;
}



#endif /* LAPLACE_STEHFEST_H_ */
