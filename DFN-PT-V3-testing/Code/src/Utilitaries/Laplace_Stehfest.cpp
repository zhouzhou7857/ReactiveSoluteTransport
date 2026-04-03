/*
 * Laplace_Stehfest.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#include "Laplace_Stehfest.h"

#include <cmath>
#include <iostream>
#include <vector>

// Stehfest multiplicator coefficients
std::vector<double> Stehfest_Calcul_Vi( int N)
{
	if(N%2!=0) std::cout << "N not even in Stehfest_Calcul_Vi\n";
	std::vector<double> V(N+1), G(N+1), H(N+1);

	// G[i] = factoriel(i)
	G[0]=1;
	for(int i=1;i<=N;i++) G[i]=G[i-1]*i;
	int Nh=N/2;
	H[1]=2/G[Nh-1];
	for(int i=2;i<=Nh;i++) H[i]=pow((double)i,(double)Nh)*G[2*i]/(G[Nh-i]*G[i]*G[i-1]);
	double sn=pow(-1.,Nh+1.);
	for(int i=1;i<=N;i++) {
		V[i]=0.;
		int lim=std::min(i,Nh);
		for(int k=(i+1)/2;k<=lim;k++)
			V[i]+=H[k]/(G[i-k]*G[2*k-i]);
		V[i]=sn*V[i];
		sn*=-1;
	}
	return V;
}

// Sampling computation
std::vector<double> Stehfest_Sampling( int N, double T) {
	std::vector<double> p(N+1);
	for( int i=1; i<=N; i++) p[i]=log(2.)/T*i;
	return p;
}

// Inverse Laplace transform by Stehfest algorithm
double Stehfet_Inverse_Laplace( double T, std::vector<double> & fp, int N, std::vector<double> & Vi) {
  double x = 0.;
  for(int i=1; i<=N; i++)
	  x+=Vi[i]*fp[i];
  x*=log(2.)/T;
  return x;
}
