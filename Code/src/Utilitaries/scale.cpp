/*
 * scale.cpp
 *
 *  Created on: Jul 6, 2012
 *      Author: delphine
 */

#include "scale.h"
#include "Constantes.h"

void scale::define( double min_, double max_, int n_, int lin_log_){
	min			= min_;
	max			= max_;
	n			= n_;
	lin_log		= lin_log_;
	if(lin_log==SCALE_LOG) {
		if(min==0) min = 1;
		d = pow(max/min,1./(n-1));
	} else
		d = (max-min)/(n-1);
}

void scale::define(std::vector<double> &scale_vect){
	n = (int) scale_vect.size();
	if(n<2)
		return;
	min = scale_vect[0];
	max = scale_vect[n-1];
	double d0 = scale_vect[1]-scale_vect[0];
	if(max-min-(n-1)*d0 < 1e-3)
		lin_log = SCALE_LIN;
	else
		lin_log = SCALE_LOG;
	if(lin_log==SCALE_LOG)
		d = pow(max/min,1./(n-1));
	else
		d = (max-min)/(n-1);
}

std::vector<double> scale::scale_vector(){
	std::vector<double> v(n);
	for(int i=0; i<n; i++)
	v[i]=scale_x_int( i);
	return v;
}

double scale::scale_x_int( int i){
	if(lin_log==SCALE_LOG)
		return (min*pow(d,i));
	else
		return (min+d*i);
}


