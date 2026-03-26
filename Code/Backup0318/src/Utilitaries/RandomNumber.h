/*
 * RandomNumber.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef RANDOMNUMBER_H_
#define RANDOMNUMBER_H_

#include <vector>
#include "../RngStream/rngstream.h"



class RngStream_a:public RngStream{
public:
	RngStream_a(){}
	~RngStream_a(){}
	double uniform(double min, double max);
	int drawPdf(std::vector<double>);
	int drawCdf(std::vector<double>&);
};

double Uniform();

#endif /* RANDOMNUMBER_H_ */
