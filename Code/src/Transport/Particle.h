/*
 * Particle.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "../Utilitaries/Pointcpp.h"
#include <vector>


class Particle {
public:
	pointcpp<double> M;	// current position of the particle
	int mesh_index;	// mesh index of the current location of the particle
	int prev_mesh_index;	// mesh index of the previous location of the particle
	std::vector<int> mesh_history;	// mesh indices visited previously
	std::vector<int> intersection_history;	// node indices of visited intersections
	double t;	// current time
	double t_injection; // injection time for travel time diagnostics
	double t_in_fract_prev; // 2026/1/2 by Wenyu, time spent in the fracture at the previous step
	int no;	// identifier of current particle
	double L_in_fract; // advective distance already covered by the particle in the current fracture
	double t_in_fract; // total time (advection+diffusion) that the particle already spent in the current fracture (since its last entrance in the fracture)
	double reactive_concentration; // reactive concentration carried by the particle; geometry updates stop when it reaches zero
	double representative_volume; // dynamic fluid volume represented by this particle, assigned at injection from inlet flux and particle time share
public:
	Particle();
	virtual ~Particle();
};

#endif /* PARTICLE_H_ */
