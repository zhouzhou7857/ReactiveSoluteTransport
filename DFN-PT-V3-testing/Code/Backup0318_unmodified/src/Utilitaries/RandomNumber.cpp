/*
 * RandomNumber.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */


#include "RandomNumber.h"
#include <iostream>
#include <stdlib.h>

using namespace std;


/** Random Number from an uniform distribution between two values : min and max */
double RngStream_a::uniform(double min, double max){
	return min + RandU01()*(max-min);
}

int RngStream_a::drawCdf(vector<double>& probas){
	double random_value = RandU01();
	int index = -1;
	for(size_t i=0 ; i<probas.size() ; i++){
		if(probas[i] > random_value){
			index = i;
			break;
		}
	}
	if(index == -1) {
		std::cerr << "WARNING in RngStream_a::drawCdf : your CDF is badly formed (not complete)\n";
		return probas.size()-1;
	}
	return index;
}

int RngStream_a::drawPdf(vector<double> probas){
	// computes the cumulative density fonction
	for(size_t i=1 ; i<probas.size() ; i++)
		probas[i] += probas[i-1];
	// draw a random number
	return drawCdf(probas);
}

double Uniform(){
	return rand()/(double)RAND_MAX;
}

