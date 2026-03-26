/*
 * BoundaryConditions.cpp
 *
 *  Created on: May 1, 2013
 *      Author: delphineroubinet
 */

#include "BoundaryConditions.h"
#include "../Utilitaries/Constantes.h"
#include <iostream>
#include <fstream>

using namespace std;


void BoundaryConditions::print(){cout << bc_type << "," << bc_value << endl;}

BoundaryConditionsDef ReadBoundaryConditions(std::string file_name,SourceTermsDef& source_terms_def){
	BoundaryConditionsDef bc_map;
	string bc_type;double bc_value,x_coord,y_coord,source_value;
	int nb_source_terms;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		// left border
		fichier >> bc_type >> bc_value;
		bc_map[LEFT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// right border
		fichier >> bc_type >> bc_value;
		bc_map[RIGHT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// top border
		fichier >> bc_type >> bc_value;
		bc_map[TOP_BORDER]=BoundaryConditions(bc_type,bc_value);
		// bottom border
		fichier >> bc_type >> bc_value;
		bc_map[BOTTOM_BORDER]=BoundaryConditions(bc_type,bc_value);
		// source terms
		fichier >> nb_source_terms;
		for (int i=0;i<nb_source_terms;i++){
			fichier >> x_coord >> y_coord >> source_value;
			source_terms_def.push_back(SourceTerm(x_coord,y_coord,source_value));
		}
	}
	else{cout << "WARNING in ReadBoundaryConditions (BoundaryConditions.cpp): Boundary conditions parameter file not found" << endl;}
	return bc_map;
}

map<string,BoundaryConditionsDef> ReadBoundaryConditionsHydroPhys(std::string file_name,SourceTermsDef& source_terms_def){
	map<string,BoundaryConditionsDef> bc_map_hydro_elec;
	BoundaryConditionsDef bc_map_elec,bc_map_hydro;
	string bc_type;double bc_value,x_coord,y_coord,source_value;
	int nb_source_terms;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		// left border
		fichier >> bc_type >> bc_value;
		bc_map_elec[LEFT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// right border
		fichier >> bc_type >> bc_value;
		bc_map_elec[RIGHT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// top border
		fichier >> bc_type >> bc_value;
		bc_map_elec[TOP_BORDER]=BoundaryConditions(bc_type,bc_value);
		// bottom border
		fichier >> bc_type >> bc_value;
		bc_map_elec[BOTTOM_BORDER]=BoundaryConditions(bc_type,bc_value);
		// left border
		fichier >> bc_type >> bc_value;
		bc_map_hydro[LEFT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// right border
		fichier >> bc_type >> bc_value;
		bc_map_hydro[RIGHT_BORDER]=BoundaryConditions(bc_type,bc_value);
		// top border
		fichier >> bc_type >> bc_value;
		bc_map_hydro[TOP_BORDER]=BoundaryConditions(bc_type,bc_value);
		// bottom border
		fichier >> bc_type >> bc_value;
		bc_map_hydro[BOTTOM_BORDER]=BoundaryConditions(bc_type,bc_value);
	}
	else{cout << "WARNING in ReadBoundaryConditions (BoundaryConditions.cpp): Boundary conditions parameter file not found" << endl;}
	bc_map_hydro_elec["elec"]=bc_map_elec;
	bc_map_hydro_elec["hydro"]=bc_map_hydro;
	return bc_map_hydro_elec;
}

BoundaryConditionsDef DefineBoundaryConditions(std::string option){
	BoundaryConditionsDef bc_map;
	if (option==CONFIG1){
		bc_map[TOP_BORDER]=BoundaryConditions(DIRICHLET,1);
		bc_map[BOTTOM_BORDER]=BoundaryConditions(DIRICHLET,0);
		bc_map[LEFT_BORDER]=BoundaryConditions(NEUMANN,0);
		bc_map[RIGHT_BORDER]=BoundaryConditions(NEUMANN,0);
	}
	else if (option==CONFIG2){
		bc_map[TOP_BORDER]=BoundaryConditions(DIRICHLET,1);
		bc_map[BOTTOM_BORDER]=BoundaryConditions(DIRICHLET,0);
		bc_map[LEFT_BORDER]=BoundaryConditions(DIR_GRADIENT);
		bc_map[RIGHT_BORDER]=BoundaryConditions(DIR_GRADIENT);
	}
	else if (option==CONFIG3){
		bc_map[LEFT_BORDER]=BoundaryConditions(DIRICHLET,1);
		bc_map[RIGHT_BORDER]=BoundaryConditions(DIRICHLET,0);
		bc_map[TOP_BORDER]=BoundaryConditions(DIR_GRADIENT);
		bc_map[BOTTOM_BORDER]=BoundaryConditions(DIR_GRADIENT);
	}
	else{cout << "WARNING in DefineBoundaryConditions (BoundaryConditions.cpp): option not defined" << endl;}
	return bc_map;
}


void print_BC(BoundaryConditionsDFN bc_map){
	for (BoundaryConditionsDFN::iterator it=bc_map.begin();it!=bc_map.end();it++){
		cout << "Boundary conditions:" << endl;
		cout << it->first << " " << it->second.bc_type << " " << it->second.bc_value << endl;
	}
}

void print_BC(BoundaryConditionsMapCurv bc_map){
	cout << "Boundary conditions:" << endl;
	for (BoundaryConditionsMapCurv::iterator it=bc_map.begin();it!=bc_map.end();it++){
		cout << it->first << " " << it->second.bc_type << " " << it->second.bc_value << endl;
	}
}


// Boundary condition definition for electric study
map<string,map<int,double> > ReturnBoundaryConditionsElec(BoundaryConditionsDef bc_map,ublas_matrix flow_velocity_x,ublas_matrix flow_velocity_y,
		map<int,double> flow_velocity_left,map<int,double> flow_velocity_bottom){
	if (bc_map[LEFT_BORDER].bc_type!=NEUMANN||bc_map[RIGHT_BORDER].bc_type!=NEUMANN||bc_map[TOP_BORDER].bc_type!=NEUMANN||bc_map[BOTTOM_BORDER].bc_type!=NEUMANN){
		cout << "WARNINIG: boundary conditions for electric study not implemented" << endl;
	}
	map<int,double> flow_velocity_right,flow_velocity_top;
	int Nx=flow_velocity_x.size1(),Ny=flow_velocity_x.size2();
	for (int j=0;j<Ny;j++){
		flow_velocity_right[j]=flow_velocity_x(Nx-1,j);
	}
	for (int i=0;i<Nx;i++){
		flow_velocity_top[i]=flow_velocity_y(i,Ny-1);
	}
	map<string,map<int,double> > flow_velocity_borders;
	flow_velocity_borders[LEFT_BORDER]=flow_velocity_left;flow_velocity_borders[RIGHT_BORDER]=flow_velocity_right;
	flow_velocity_borders[TOP_BORDER]=flow_velocity_top;flow_velocity_borders[BOTTOM_BORDER]=flow_velocity_bottom;
	return flow_velocity_borders;
}
