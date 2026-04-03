/*
 * FlowComputation.cpp
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */

#include "FlowComputation.h"
#include "../Utilitaries/UblasStructures.h"
#include "../Utilitaries/LinearSystem.h"
#include "../Domain_Definition/HydraulicProperties.h"

using namespace std;


/**********************************/
// Functions for Flow Computation //
/**********************************/

/*
ublas_vector DFNComputationClassic(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map){
	// 1. Nodes definition
	NodeFracturesMap node_fract_map=net_mesh.ReturnNodeFracturesMap();
	int nb_nodes=node_fract_map.size();
	// 2. Linear system definition
	ublas_matrix mat_syst(nb_nodes,nb_nodes); ublas_vector vect_syst(nb_nodes);
	// 2.1. Loop on each node
	double value;
	for (NodeFracturesMap::iterator it1=node_fract_map.begin();it1!=node_fract_map.end();it1++){
		// Matrix affectation - Loop on corresponding fractures
		for (HydrauPropMap::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			value=it2->second.transmissivity/it2->second.length;
			mat_syst(it1->first,it2->first)=-value;
			mat_syst(it1->first,it1->first)+=value;
		}
		// Vector affectation
		vect_syst(it1->first)=0;
	}
	// 3. Boundary conditions affectation
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// cancel previous values
		for (int i=0;i<mat_syst.size2();i++){mat_syst(it1->first,i)=0;}
		// new value affectation
		mat_syst(it1->first,it1->first)=1;
		vect_syst(it1->first)=it1->second;
	}
	// 4. Linear system solution
	return LinearSystemSolving(mat_syst,vect_syst);
}

BoundaryConditionsDFN ReturnBoundCondDFN(Domain domain,BordersMap border_map,NodesMap nodes_map,double left_pressure,double right_pressure){
	// 1. Variables
	BoundaryConditionsDFN bc_map;
	// 2. Loop on each node of the fracture network
	for (map<int,string>::iterator it=border_map.begin();it!=border_map.end();it++){
		if (it->second==LEFT_BORDER){bc_map[it->first]=left_pressure;}
		else if (it->second==RIGHT_BORDER){bc_map[it->first]=right_pressure;}
		else{
			string border1,border2;
			double Lx=domain.domain_size_x();
			double position=domain.ReturnBorderCoordinate(nodes_map[it->first],Lx,border1,border2);
			if (border1!=LEFT_BORDER||border2!=RIGHT_BORDER){
				bc_map[it->first]=(right_pressure-left_pressure)/Lx*position+left_pressure;
			}
		}
	}
	return bc_map;
}

double VelocityComputation(double aperture,double h1,double h2,double length){
	double transmissivity=ReturnTransmissivity(aperture);
	return transmissivity/aperture*(h1-h2)/length;
}

// Compute flow and return the fracture network composed of fracture with flow
NetworkMeshes FlowComputation(NetworkMeshes net_mesh){
	NetworkMeshes flow_net_mesh=net_mesh;
	// 1. Boundary condition definition
	BoundaryConditionsDFN bc_map=ReturnBoundCondDFN(net_mesh.domain,net_mesh.border_map,net_mesh.nodes_map,net_mesh.left_pressure,net_mesh.right_pressure);
	// 2. Pressure computation
	ublas_vector hydraulic_head=DFNComputationClassic(net_mesh,bc_map);
	// 3. Flow computation
	double velocity,max_velocity=0;FractureMesh fract_mesh;
	vector<FractureMesh> new_meshes;
	for (int i=0;i<net_mesh.meshes.size();i++){
		velocity=VelocityComputation(net_mesh.meshes[i].aperture,hydraulic_head(net_mesh.meshes[i].p_ori.index),hydraulic_head(net_mesh.meshes[i].p_tar.index),net_mesh.meshes[i].ReturnLength());
		if (fabs(velocity)>EPSILON){
			if (velocity>0){
				fract_mesh=net_mesh.meshes[i];
				fract_mesh.velocity=velocity;
			}
			else if (velocity<0){
				fract_mesh=net_mesh.meshes[i];
				fract_mesh.p_ori=net_mesh.meshes[i].p_tar;
				fract_mesh.p_tar=net_mesh.meshes[i].p_ori;
				fract_mesh.velocity=-velocity;
			}
			new_meshes.push_back(fract_mesh);
			if (velocity>max_velocity){max_velocity=velocity;}
		}
	}
	flow_net_mesh.meshes=new_meshes;
	cout << max_velocity << endl;
	return flow_net_mesh;
}
*/
