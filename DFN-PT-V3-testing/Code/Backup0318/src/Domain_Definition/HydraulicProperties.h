/*
 * HydraulicProperties.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef HYDRAULIC_PROPERTIES_H_
#define HYDRAULIC_PROPERTIES_H_

#include "../Domain_Definition/FractureMesh.h"

class HydraulicProperties{
public:
	double length;
	double aperture;
	double transmissivity;
public:
	HydraulicProperties();
	virtual ~HydraulicProperties();
	HydraulicProperties(FractureMesh);
	double ReturnTransmissivity();
};

double ReturnTransmissivity(double);

#endif /* HYDRAULIC_PROPERTIES_H_ */
