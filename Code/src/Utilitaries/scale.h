/*
 * scale.h
 *
 *  Created on: Jul 6, 2012
 *      Author: delphine
 */

#ifndef SCALE_H_
#define SCALE_H_

#include <vector>

class scale {
private:
	double	min;	// Minimal class value
	double	max;	// Maximal class value
	double	d;		// Increment of class
	int		n;		// Number of classes
	int		lin_log;// lin: 0, log=1
	std::vector<double> scale_vect_;
public:
	scale(double min_, double max_, int n_, int lin_log_){this->define(min_,max_,n_,lin_log_);}
	scale(){ min=0.; max=0.; d = 0.; n=0; lin_log= 0;};
	~scale(){};
	scale(const scale & o){ min = o.min; max = o.max; n = o.n; d = o.d; lin_log = o.lin_log;}
	scale& operator=(const scale & o){ min = o.min; max = o.max; n = o.n; d = o.d; lin_log = o.lin_log; return *this; }

	void define( double min_, double max_, int n_, int lin_log_);
	void define(std::vector<double> &scale_vect);
	int scale_int_x( double x);
	double scale_x_int( int i);
	std::vector<double> scale_vector();
};

#endif /* SCALE_H_ */
