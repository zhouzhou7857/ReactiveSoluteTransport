/*
 * Transport.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "../Input_Output/Parameters.h"
#include "../Domain_Definition/Domain.h"
#include "../Domain_Definition/NetworkMeshes.h"
#include "../Transport/Particle.h"
#include "../Utilitaries/RandomNumber.h"
#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Segment.h"
#include "../Utilitaries/Point_Cgal.h"
#include "../Utilitaries/Structures.h"

#include <boost/math/special_functions/erf.hpp>

class PhysicsParam{
public:
	double Dm;	// matrix diffusion coefficient
	double porosity;	// matrix porosity
public:
	PhysicsParam(){};
	virtual ~PhysicsParam(){};
	PhysicsParam(double Dm_,double porosity_){Dm = Dm_; porosity = porosity_;};
};

class NumericalParam{
public:
	int nb_part;	// number of particles for transport simulation
	double proba_transfer;	// determine segment discretization for particle transport
	int simu_option;
	double t_max;
	double t_injection;
	double output_interval;
	double reaction_dt;
public:
	NumericalParam(){};
	virtual ~NumericalParam(){};
	NumericalParam(int nb_part_,double proba_transfer_,int simu_option_,double t_max_,double t_injection_,double output_interval_,double reaction_dt_){
		nb_part = nb_part_;
		proba_transfer = proba_transfer_;
		simu_option = simu_option_;
		t_max = t_max_;
		t_injection = t_injection_;
		output_interval = output_interval_;
		reaction_dt = reaction_dt_;
	};
};

struct ParticleSnapshot{
	int no;
	int mesh_index;
	int mesh_index_original;
	double t_injection;
	pointcpp<double> M;
};

class Transport{
public:
	RngStream_a rng_tracker;		///< Random Number Generator for diffusion and dispersion
	NetworkMeshes net_mesh;	// description of the fracture network
	NetworkMeshes net_mesh_modified;  // added by DR on 2025/12/11: network with aperture meshes modified due to chemical reactions
	NetworkMeshes net_mesh_before_backbone; // snapshot before backbone rebuild for output comparison
	Domain domain;	// geometric description of the domain
	PhysicsParam phys_param;
	NumericalParam num_param;
	Parameters param_full;
	Distribution Transfer_Time_Distribution;	///< Transfer time distributions
	std::map<int,pointcpp<double> > Initial_Position_Particles;
	std::map<double,std::vector<ParticleSnapshot> > Position_Time_Particles;
	MapSegmentTimes Segment_Particle_Time; // For each segment, we store the time at which each particle leaves (or reaches another segment)
	MapSegmentVelocities Segment_Velocities; // For each segment, we store the flow velocities
	Declare_CDF_Map Delta_CDF_Map; // pair(velocity,abs(sigma_end-sigma_beg),delta_cdf;
	Projection_Map_Global proj_map_global;	// store the projected points associated with the pair of points defined as keys of the map
public:
	Transport(){};
	virtual ~Transport(){};
	Transport(RngStream_a,NetworkMeshes,Domain,Parameters);
	std::vector<Particle> Particles_Injection(int option_injection=0);
	bool infinite_matrix_displacement(Particle&);
	bool finite_matrix_displacement(Particle&);
	bool infinite_matrix_displacement_step(Particle&,double,std::map<int,double>&);
	bool finite_matrix_displacement_step(Particle&,double,std::map<int,double>&);
	bool particle_displacement(Particle&);
	bool particle_displacement_step(Particle&,double,std::map<int,double>&);
	bool Particles_Transport(std::map<int,double>&,int option_injection=0);
	void RecordPositions(double,const std::vector<Particle>&);
	void WritePositionSnapshotsCSV(const std::string&) const;
	void PrintTravelTimeSamplesAtOutputIntervals(int max_samples=5) const;
	double Get_Total_Time_From_Advec_Time(double,double,double,double);
	//int Mesh_Neighbour_Choice_And_Test(pointcpp<double>,FractureMesh,RngStream_a&);
	bool Mesh_Neighbour_Choice_And_Test(pointcpp<double>,FractureMesh,RngStream_a&,int&,bool&);
	bool select_neighbour_slr(std::vector<int>&,std::vector<double>&,FractureMesh,CgalPoint2D);

	// functions for orthogonal projection between fractures
	bool Projection_Periodic_Boundaries(CgalPoint2D,Segment2D,Segment2D,CgalPoint2D,CgalPoint2D,CgalPoint2D&,int&,double&);
	bool Closest_Orthogonal_Projection(CgalPoint2D,Segment2D,CgalPoint2D&,int&,double&);
	//Projection_Map Orthogonal_Projection_M(CgalPoint2D,CgalPoint2D);
	bool Orthogonal_Projection_M(CgalPoint2D,CgalPoint2D,Projection_Map&);
	double ReturnPureDiffusionTime(double);
	bool ReflectionOnBorder(pointcpp<double>&,pointcpp<double>&,int &,double &);
	void M_project_in_system(CgalPoint2D&,const double&);
	void M_project_in_system(pointcpp<double> &,const double &);

	// functions for transferring to neighboring fractures
	bool Transfer_Probability_and_Transfer_Time(Projection_Map,const double,int&,double&,pointcpp<double>&);
	bool Transfer_Probability_Computation(const std::pair<double,double> &,const double,double &);
	bool Transfer_Time_Computation(const std::pair<double,double>&,const double,double&);

	bool exist_in_map(pointcpp<double>,pointcpp<double>,Projection_Map &);
};

// functions for transferring probability
double Transfer_Probability_FPTD(const double,const double,const double);
bool Transfer_Probability_Feller(const std::pair<double,double> &,const double,const double,double&);
bool Feller_Analytical_Solution(double&,const double,const std::pair<double,double>&,const double);
bool Feller_Analytical_Solution_Single(double&,const double,const double,const double,const double);
double Transfer_Time_FPTD(const double,const double,const double);
bool Transfer_Time_Feller(const std::pair<double,double>&,const double,const double,RngStream_a&,double&,Distribution&);
bool Transfer_Time_Distribution_Computation(Distribution &,const std::pair<double,double> &,const double,const double,RngStream_a&);
bool Feller_Maximal_Time_Determination(double &,const double,const std::pair<double,double> &,const double);
bool Get_Time_Feller(std::pair<std::vector<double>,std::vector<double> >&,const double,double&);
bool Barrier_Choice(const std::map<int,double>&,const double,int&);
bool Barrier_Choice(std::pair<double,double>, const double, std::string &);




#endif /* TRANSPORT_H_ */
