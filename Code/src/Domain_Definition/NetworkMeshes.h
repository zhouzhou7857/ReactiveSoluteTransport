/*
 * NetworkMeshes.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef NETWORKMESHES_H_
#define NETWORKMESHES_H_

#include <set>
#include "../Input_Output/Parameters.h"
#include "FractureMesh.h"
#include "Domain.h"
#include "HydraulicProperties.h"
#include "../Utilitaries/UblasStructures.h"
#include "../Chemistry/Chemistry.h"

typedef std::map<int,double> SourceTermDFN;	//<node index,source term>

typedef std::map<int,std::string> BordersMap; // <node index,border name>
typedef std::map<int,CgalPoint2D> NodesMap;

typedef std::map<int,HydraulicProperties> HydrauPropMap;
typedef std::map<int,HydrauPropMap> NodeFracturesMap;

typedef std::map<int,std::pair<FluxPoint2D,FluxPoint2D> > FracturesMap;


class NetworkMeshes {
public:
	Domain domain;
	std::vector<FractureMesh> meshes;	// list of meshes composing the fracture network
	int cpt_inter;	// fracture nodes (intersections and extremities) numbering
	int cpt_mesh;	// fracture mesh numbering
	int cpt_fract; // fracture counter (index of the last fracture, nb fract=cpt_fract+1 as index starts at 0)
	BordersMap border_map;	// map<node_index,domain_border>
	NodesMap nodes_map;	// map<node_index,node_coordinate>
	std::set<int> inter_list;
	double max_fract_spacing;
	// map of fractures with the index of each fracture and their extremities (origin,target) with a positive flow velocity going from the origin to the target of the fractures
	std::map<int,std::pair<FluxPoint2D,FluxPoint2D> > fractures;
	std::set<int> fract_indices; // set of the indices of the fractures which is used to define the fracture map
public:
	NetworkMeshes();
	virtual ~NetworkMeshes();
	NetworkMeshes(std::string,std::string,Domain);
	//NetworkMeshes(Domain,double,double,double,double,double,double,int);
	NetworkMeshes(Domain,double,double,double,double,int);
	NetworkMeshes(double,double,Parameters,Domain);
	FractureMesh return_mesh(int);
	FractureMesh return_mesh(int) const;
	std::vector<FractureMesh> return_mesh_from_node_ori(int);
	std::vector<int> Nb_passed;   // Nb for each mesh: number of particles passed through mesh i
	std::vector<double> last_t_update;
	FractureMesh operator[](int index){return meshes[index];}
	int ReturnIndex(CgalPoint2D);
	void print_DFN();
	void print_fractures();
	void print_DFN_in_file(Parameters);
	void print_DFN_in_file(std::string,int);
        void print_DFN_in_file_initial(std::string,int);
        void print_DFN_in_file_initial(Parameters);
	void print_DFN_in_file_with_aperture_delta(Parameters,const NetworkMeshes&);
	NodeFracturesMap ReturnNodeFracturesMap();
	std::set<int> return_connected_nodes();
	NetworkMeshes remove_fractures_from_velocity(Parameters,double,std::string);
	NetworkMeshes return_network_from_nodes(std::set<int>);
	std::set<int> return_nodes_from_velocity(Parameters,double,std::string option="neumann");
	NetworkMeshes return_backbone(Parameters,double,std::string option="neumann");
	void EvaluateFlowVelocities(ublas_vector);
	std::pair<double,double> return_min_max_vel();
	double return_ave_vel();
	bool define_fracture_map();
	double ReturnVolumetricFractureDensity();
	void print_DFN_in_file_aperture(std::string,int);
	// Legacy empirical aperture update interface.
	// 在 DFN-PT-V3 当前主路径中并未使用。
	void ChangeAperture(int, double, double, double);
	// Legacy empirical ratio-based aperture update interface.
	// 在 DFN-PT-V3 当前主路径中并未使用。
	void ChangeApertureByRatio(int, double);
};

void NodeInsertionMap(int,int,FractureMesh,NodeFracturesMap&);
void add_neigh_nodes(std::set<int>&,NodeFracturesMap,int);
void NetworkMeshesFromFile(std::string,std::string,NetworkMeshes&);
void NetworkMeshesGenerationSierpinski(std::string,std::string,NetworkMeshes&);
void NetworkMeshesGenerationRealistic(std::string,std::string,NetworkMeshes&);
void NetworkMeshesGenerationRealistic2(std::string,std::string,NetworkMeshes&);
void NetworkMeshesGenerationRealistic3(std::string,std::string,NetworkMeshes&);
void DivideDomain(pointcpp<double>,pointcpp<double>,int,double,NetworkMeshes&,double);
std::vector<Domain> DivideSierpinskiDomain(pointcpp<double>,int,double,double,double,NetworkMeshes&,double);
void AddFracture(FractureMesh,NetworkMeshes&);
void Border_Definition(FluxPoint2D,NetworkMeshes&);
void ComputeIntersections(NetworkMeshes&);
NetworkMeshes Vertical_Translation(NetworkMeshes,double);
void FractureMeshNumbering(NetworkMeshes&);
void UpdateCptFract(NetworkMeshes&);
//void VerticalUpdateFracture(NetworkMeshes&,std::vector<FractureMesh>&,int,double,double);


/***** FUNCTIONS TO COMPUTE THE FLOW VELOCITY *******/
double VelocityComputation(double,double,double,double);
NetworkMeshes FlowComputation(NetworkMeshes,BoundaryConditionsDef);


#endif /* NETWORKMESHES_H_ */
