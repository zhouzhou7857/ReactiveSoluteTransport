/*
 * Parameters.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

//enum{INFINITE_MATRIX,FINITE_MATRIX};

#include <string>
#include <set>

#include "../Utilitaries/Constantes.h"
#include "../Input_Output/BoundaryConditions.h"

typedef std::map<std::pair<pointcpp<double>,pointcpp<double> >,std::set<double> > MapSegmentTimes;
typedef std::map<std::pair<pointcpp<double>,pointcpp<double> >,double > MapSegmentVelocities;

class Parameters{
public:
	double Lx;	// domain size in the longitudinal direction
	double Ly;	// domain size in the transversal direction
	double Dm;	// matrix diffusion coefficient
	double porosity;	// matrix porosity
	long int nb_part;	// number of particles for transport simulation
	double proba_transfer;	// determine segment discretization for particle transport
	int simu_option;	// 0/1 = INFINITE_MATRIX/FINITE_MATRIX
	int Nt;	// number of time steps for result post-processing
	double t_min;	// minimum time for result post-processing
	double t_max;	// maximum time for result post-processing
	double t_injection;	// duration of particle injection; Modified by Wenyu on 2026/1/6
	double output_interval;	// time interval for position outputs; Modified by Wenyu on 2026/1/8
	double reaction_dt;	// fixed reaction update time step used for transport/aperture updates
	double mineral_volume_reference_water_volume;	// PHREEQC reference water volume Vref
	double fracture_out_of_plane_thickness;	// thickness used in volume-to-aperture conversion
	//int num_simu_start;	// number of the starting simulation
	std::string code_path;
	std::string file_names_path;
	std::string file_name_domain,file_name_DFN,file_name_simu,file_name_chemistry;	// parameter files
	BoundaryConditionsDef bc_map_def;	// boundary conditions
	//int seed;	// seed used for the generating the fracture networks (same seed for each pair of (C,D))
	// DFN parameters
	std::string generation_option_DFN;
	double density_param,exponent_param,fract_aperture,r_min,coeff_theta1,coeff_theta2,b_min,b_max,mean_lnb,RSD_lnb;
	int seed_simu;
public:
	Parameters(const std::string& file_names_override = "");
	~Parameters(){};
	void read_param();
};

class InvParam{
public:
	std::vector<double> density_param;	// fracture density
	std::vector<double> exponent_param;	// fractal dimension
	std::vector<double> coeff_theta1;		// first preferential angle
	std::vector<double> coeff_theta2;		// second preferential angle
public:
	InvParam();
	~InvParam(){};
};

#endif /* PARAMETERS_H_ */
