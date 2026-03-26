/*
 * DFNComputation.cpp
 *		HEAT
 *  	Created on: August 31, 2015
 *      Author: viktoria and delphine
 */


#include "DFNComputation.h"
#include "../Utilitaries/LinearSystem.h"

using namespace std;

BoundaryConditionsDFN ReturnBoundCondDFN(Domain domain,BoundaryConditionsDef bc_map_classic,BordersMap border_map,NodesMap nodes_map){
	// 1. Variables
	BoundaryConditionsDFN bc_map;int index;string bc_type,border1,border2;double position,border_dist,bc_value;
	// 2. Loop on each node of the fracture network
	for (map<int,string>::iterator it=border_map.begin();it!=border_map.end();it++){
		index=it->first;bc_type=bc_map_classic[it->second].bc_type;
		// 2.1. Classical affectation (boundary condition value is known)
		if (bc_type==DIRICHLET||bc_type==NEUMANN){
			bc_map[index]=bc_map_classic[it->second];
		}
		// 2.2. Specific affectation (boundary condition value has to be defined)
		else if (bc_type==DIR_GRADIENT){
			// Position of the node on the border
			position=domain.ReturnBorderCoordinate(nodes_map[index],border_dist,border1,border2);
			// Boundary condition value at this position
			if (bc_map_classic[border1].bc_type!=DIRICHLET||bc_map_classic[border2].bc_type!=DIRICHLET){
				cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): case not implemented" << endl;
			}
			else{
				bc_value=(bc_map_classic[border2].bc_value-bc_map_classic[border1].bc_value)/border_dist*position+bc_map_classic[border1].bc_value;
				bc_map[index]=BoundaryConditions(DIRICHLET,bc_value);
			}
		}
		else if (bc_type==MIXED_BC){cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): boundary condition type not implemented" << endl;}
		else{cout << "WARNING in ReturnBoundCondDFN (DFNComputation.cpp): boundary condition type not implemented" << endl;}
	}
	return bc_map;
}

ublas_vector DFNComputationClassic(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map){
	SourceTermDFN source_term_DFN;
	return DFNComputationClassic(net_mesh,bc_map,source_term_DFN);
}

ublas_vector DFNComputationClassic(NetworkMeshes net_mesh,BoundaryConditionsDFN bc_map,SourceTermDFN source_term_DFN,string option,set<int> connected_nodes){



		//ublas_vector DFN_potential){
	// 1. Nodes definition
	NodeFracturesMap node_fract_map=net_mesh.ReturnNodeFracturesMap();
	int nb_nodes=node_fract_map.size();
	// 2. Linear system definition
	ublas_matrix mat_syst(nb_nodes,nb_nodes); ublas_vector vect_syst(nb_nodes);

	// 2.1. Loop on each node
	double value;
	for (NodeFracturesMap::iterator it1=node_fract_map.begin();it1!=node_fract_map.end();it1++){
		// Matrix affectation - Loop on corresponding fractures  // it2 is each node connected to the first node it1 we run through it2
		for (HydrauPropMap::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			value=it2->second.ReturnTransmissivity()/it2->second.length; // defines the coefficient in the A matrix
			mat_syst(it1->first,it2->first)=-value;  // for p1 it is plus this value
			mat_syst(it1->first,it1->first)+=value;  // for p2 it is minus this value
		}
		// Vector affectation  vector b all the flow rate is equal zero so we sett it here
		vect_syst(it1->first)=0;
	}

	// 3. Source term affectation
	for (SourceTermDFN::iterator it=source_term_DFN.begin();it!=source_term_DFN.end();it++){
		vect_syst(it->first)+=it->second;
	}

	// 4. Boundary conditions affectation
	HydrauPropMap elec_prop;
	for (BoundaryConditionsDFN::iterator it1=bc_map.begin();it1!=bc_map.end();it1++){
		// electric potential affectation
		if (it1->second.bc_type==DIRICHLET){
			// cancel previous values
			elec_prop=node_fract_map[it1->first];
			for (HydrauPropMap::iterator it2=elec_prop.begin();it2!=elec_prop.end();it2++){mat_syst(it1->first,it2->first)=0;}
			// new value affectation
			mat_syst(it1->first,it1->first)=1;
			vect_syst(it1->first)=it1->second.bc_value;
		}
		// electric current affectation
		else if (it1->second.bc_type==NEUMANN){
			// new value affectation
			vect_syst(it1->first)=it1->second.bc_value;
		}
		else{cout << "WARNING in DFNComputation.cpp: boundary condition not implemented" << endl;}
	}
	// 5. Linear system solution  --- solve x aka p
	return LinearSystemSolving(mat_syst,vect_syst);
}



ublas_vector ComputeFlowVelocities(NetworkMeshes& net_mesh,Parameters param,Domain domain,string option){
	// Definition of DFN boundary conditions
	if (option=="gradient"){
		param.bc_map_def[TOP_BORDER]=BoundaryConditions(DIR_GRADIENT,0.);
		param.bc_map_def[BOTTOM_BORDER]=BoundaryConditions(DIR_GRADIENT,0.);
	}
	else if (option=="neumann"){
		param.bc_map_def[TOP_BORDER]=BoundaryConditions(NEUMANN,0.);
		param.bc_map_def[BOTTOM_BORDER]=BoundaryConditions(NEUMANN,0.);
	}
	else{cout << "WARNING in DFNComputation.cpp:option not implemented" << endl;}
	BoundaryConditionsDFN bc_map_dfn=ReturnBoundCondDFN(domain,param.bc_map_def,net_mesh.border_map,net_mesh.nodes_map);
	// hydraulic head computation
	SourceTermDFN source_term_DFN;
	ublas_vector DFN_potential=DFNComputationClassic(net_mesh,bc_map_dfn,source_term_DFN,"velocity");
	// flow velocity computation
	net_mesh.EvaluateFlowVelocities(DFN_potential);
	return DFN_potential;
}
