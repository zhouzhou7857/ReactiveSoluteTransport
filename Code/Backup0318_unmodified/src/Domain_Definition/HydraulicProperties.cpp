/*
 * HydraulicProperties.cpp
 *		HEAT
 *  	Created on: August 31, 2015
 *      Author: viktoria and delphine
 */

#include "HydraulicProperties.h"

using namespace std;

/************************************/
// Functions for Hydraulic Properties //
/************************************/
HydraulicProperties::HydraulicProperties(){}
HydraulicProperties::~HydraulicProperties(){}
HydraulicProperties::HydraulicProperties(FractureMesh fract_mesh){
	length=fract_mesh.ReturnLength();
	aperture=fract_mesh.aperture;
	transmissivity=ReturnTransmissivity();
}

// T=rho*g*b^3/(12*mu) is approximated by T=rho*g*b^3/(10*mu) with g approximated by 10
// (to easily enforced known value of the velocity)
double HydraulicProperties::ReturnTransmissivity(){
	return RHO*G*std::pow(aperture,3)/(12*MU);
}


double ReturnTransmissivity(double aperture){
	return RHO*G*std::pow(aperture,3)/(12*MU);
}
