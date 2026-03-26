/*
 * Results.h
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */

#ifndef RESULTS_H_
#define RESULTS_H_

#include "../Input_Output/Parameters.h"
#include "../Utilitaries/Structures.h"
#include <map>
#include <string>


class Results{
public:
	int Nt;	// number of time steps for results post-processing
	double t_min,t_max;	// minimum and maximum times for results post-processing
	std::map<int,double> arrival_times;	//<particle identifier,arrival time>
	std::map<double,int> cum_dist_times;	//<arrival time,number of particles arrived before this time>
	std::map<double,double> pdf_part;	//<arrival time,number of particles arrived before this time>
public:
	Results(){Nt=0;t_min=0.0;t_max=0.0;};
	virtual ~Results(){};
	Results(std::map<int,double>,Parameters);
	void post_processing();
	void writing(std::string);
	void writing(std::string,int);
};

pointcpp<double> ReadCoordinates(std::string);

#endif /* RESULTS_H_ */
