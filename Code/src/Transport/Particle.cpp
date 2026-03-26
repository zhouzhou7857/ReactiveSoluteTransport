/*
 * Particle.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Particle.h"

Particle::Particle(){
	M=pointcpp<double>(-1,-1,-1);
	mesh_index=-1;
	prev_mesh_index=-1;
	mesh_history.clear();
	intersection_history.clear();
	t=-1;
	t_injection=-1;
	no=-1;
	L_in_fract=0;
	t_in_fract=-1;
	t_in_fract_prev = 0.0;// 2026/1/2 by Wenyu, time spent in the fracture at the previous step
}
Particle::~Particle(){}
