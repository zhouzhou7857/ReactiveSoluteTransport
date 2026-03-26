/*
 * Parameters.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Parameters.h"
#include <fstream>
#include <iostream>

using namespace std;

/********************** PARAMETERS **************************/
Parameters::Parameters(){
	code_path="./../..";
	//code_path="/Users/delphineroubinet/Documents/Projets/FracTherm/CodeHeatTransport";
	string file_name = code_path+"/Input/File_names.txt";
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> file_name_domain;
		fichier >> file_name_simu;
		fichier >> file_name_DFN;
	}
	else{cout << "Parameters files not found" << endl;}
	read_param();
}


// Parameter reading and affectation
void Parameters::read_param(){
	// 1. Read the domain parameters and build the corresponding boundary conditions
	// Read the parameters
	string file_name=code_path+"/Input/Domain_files/"+file_name_domain;
	ifstream fichier1(file_name.c_str());
	double left_head=0.0,right_head=0.0;
	if (fichier1.is_open()){
		fichier1 >> Lx >> Ly;
		fichier1 >> Dm >> porosity;
		fichier1 >> left_head >> right_head;
	}
	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file1" << endl;
	}
	fichier1.close();
	// Build the boundary conditions
	bc_map_def[LEFT_BORDER]=BoundaryConditions(DIRICHLET,left_head);
	bc_map_def[RIGHT_BORDER]=BoundaryConditions(DIRICHLET,right_head);
	bc_map_def[BOTTOM_BORDER]=BoundaryConditions(NEUMANN,0);
	bc_map_def[TOP_BORDER]=BoundaryConditions(NEUMANN,0);

	// 2. Read the simulation parameters and define the matrix option
	// Read the parameters
	file_name=code_path+"/Input/Simulation_files/"+file_name_simu;
	ifstream fichier2(file_name.c_str());
	if (fichier2.is_open()){
		fichier2 >> nb_part;
		fichier2 >> proba_transfer;
		fichier2 >> simu_option;
		fichier2 >> t_min;
		fichier2 >> t_max;
		fichier2 >> Nt;
		//fichier2 >> num_simu_start;
		fichier2 >> seed_simu;
		if (!(fichier2 >> t_injection)){
			t_injection = t_max;
		}
		if (!(fichier2 >> output_interval)){
			output_interval = t_injection;
		}
		if (!(fichier2 >> reaction_dt)){
			if (nb_part>1){
				reaction_dt = t_injection/(double)(nb_part-1);
			}
			else{
				reaction_dt = t_injection;
			}
		}
	}
	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file2" << endl;
	}
	fichier2.close();

	// 3. Read the parameters related to the fracture network
	file_name=code_path+"/Input/DFN_files/"+file_name_DFN;
        ifstream fichier3(file_name.c_str());
        if (fichier3.is_open()){
                fichier3 >> generation_option_DFN;
		/*if (generation_option_DFN=="generation_realistic2"){
			fichier3 >>  density_param >> exponent_param >> fract_aperture >> r_min >> coeff_theta1 >> coeff_theta2;
		}
		else if (generation_option_DFN=="generation_realistic3"){
			fichier3 >>  density_param >> exponent_param >> b_min >> b_max >> mean_lnb >> RSD_lnb >> r_min;
                }
		else{cout << "WARNING in Parameters (Parameters.cpp) : option_generation_DFN not implemented" << endl;}*/
        }
        else{cout << "WARNING in Parameters (Parameters.cpp) : unknown file2" << endl;}
        fichier2.close();
}

/********************** INV_PARAM **************************/
InvParam::InvParam(){
	//string file_name = "./../../Input/C_d.txt";
	string file_name = "./../../Input/p_a.txt";
	cout << "Read the input inversion parameters file: " << file_name << endl;
	double value;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		while (!fichier.eof()){
			fichier >> value; 
			if (!fichier.eof()){
				density_param.push_back(value);
				fichier >> value; exponent_param.push_back(value);
				/*fichier >> value; coeff_theta1.push_back(value);
				fichier >> value; coeff_theta2.push_back(value);*/
			}
		}
	}
	else{cout << "Parameters files not found" << endl;}
}
