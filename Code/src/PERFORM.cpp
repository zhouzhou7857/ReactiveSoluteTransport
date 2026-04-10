//============================================================================
// Name        : PERFORM.cpp

// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : C++ code for solute transport in heterogeneous fractured porous media
//============================================================================

#include "Input_Output/Parameters.h"
#include "Domain_Definition/Domain.h"
#include "Domain_Definition/NetworkMeshes.h"
#include "Transport/Transport.h"
#include "Input_Output/Results.h"
#include "Chemistry/Chemistry.h"
#include "Utilitaries/Structures.h"
#include "Visualisation/DFNVisu.h"
#include "Visualisation/DisplayResults.h"
#include "Flow/FlowComputation.h"
#include <iostream>
#include <chrono>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {

	clock_t t_begin_total=clock(),t_begin_simu;

	int option=1; // 0/1: fracture network generation/transport simulation
	int option_injection=2; // 0: extended-uniform / 1: localized injection / 2: extended-flux-weighted
	bool backbone=true;

	// 1. Read parameters and define domain
	// 1.1. Standard definition
	Parameters param;
	ConfigureChemistryParameters(param.chemistry_initial_reactive_concentration,
			param.chemistry_reactive_concentration_decay,
			param.chemistry_reactive_to_mineral_stoich,
			param.chemistry_mineral_molar_volume,
			param.chemistry_fracture_out_of_plane_thickness);
	Domain domain(param.Lx,param.Ly);

	// 1.2. Parameters related to the fracture network
	cout << "0. Parameters" << endl;
	cout << "fract_aperture generation = " << param.generation_option_DFN << endl;
	cout << "nb part = " << param.nb_part << endl;	
	cout << "p_lim = " << param.proba_transfer << endl;
	cout << "option matrix = " << param.simu_option << endl;
	cout << "DFN file = " << param.file_name_DFN << endl;
	cout << "Chemistry file = " << (param.file_name_chemistry.empty() ? std::string("(default)") : param.file_name_chemistry) << endl;
	cout << "Chemistry parameters: C0=" << INITIAL_REACTIVE_CONCENTRATION
	     << " decay=" << REACTIVE_CONCENTRATION_DECAY
	     << " stoich=" << REACTIVE_TO_MINERAL_STOICH
	     << " molar_volume=" << MINERAL_MOLAR_VOLUME
	     << " thickness=" << FRACTURE_OUT_OF_PLANE_THICKNESS << endl;
	cout << "Domain bounds: min(" << domain.min_pt.i << "," << domain.min_pt.j
	     << "), max(" << domain.max_pt.i << "," << domain.max_pt.j << ")" << endl;

	t_begin_simu=clock();
	// 2. Fracture network generation and backbone definition
	cout << "1. Fracture network generation" << endl;
	auto connectivity_start = std::chrono::steady_clock::now();
	NetworkMeshes init_net_mesh;
	if (param.generation_option_DFN=="generation_realistic3"
	 || param.generation_option_DFN=="generation_realistic4"){
		cout << "DFN generation path: constructor-based " << param.generation_option_DFN << endl;
		init_net_mesh=NetworkMeshes(param.density_param,param.exponent_param,param,domain);
	}
	else{
		init_net_mesh=NetworkMeshes(param.code_path,param.file_name_DFN,domain);
	}
	NetworkMeshes net_mesh=init_net_mesh;
	int mesh_count_before_backbone = init_net_mesh.meshes.size();

	auto count_border_meshes = [&](const NetworkMeshes& mesh_net){
		int left_count = 0;
		int right_count = 0;
		for (size_t i=0;i<mesh_net.meshes.size();i++){
			pointcpp<double> p_ori(mesh_net.meshes[i].p_ori.p.x(),mesh_net.meshes[i].p_ori.p.y());
			pointcpp<double> p_tar(mesh_net.meshes[i].p_tar.p.x(),mesh_net.meshes[i].p_tar.p.y());
			if (domain.on_input_limit(p_ori) || domain.on_input_limit(p_tar)){left_count++;}
			if (domain.on_output_limit(p_ori) || domain.on_output_limit(p_tar)){right_count++;}
		}
		return make_pair(left_count,right_count);
	};
	auto count_in_domain_meshes = [&](const NetworkMeshes& mesh_net){
		int in_domain = 0;
		for (size_t i=0;i<mesh_net.meshes.size();i++){
			if (domain.IsInDomain(mesh_net.meshes[i].p_ori.p) && domain.IsInDomain(mesh_net.meshes[i].p_tar.p)){
				in_domain++;
			}
		}
		return in_domain;
	};
	auto report_dfn_extents = [&](const NetworkMeshes& mesh_net){
		if (mesh_net.meshes.empty()){
			cout << "DFN extents: meshes empty" << endl;
			return;
		}
		double min_x = mesh_net.meshes[0].p_ori.p.x();
		double max_x = min_x;
		double min_y = mesh_net.meshes[0].p_ori.p.y();
		double max_y = min_y;
		for (size_t i=0;i<mesh_net.meshes.size();i++){
			double x1 = mesh_net.meshes[i].p_ori.p.x();
			double y1 = mesh_net.meshes[i].p_ori.p.y();
			double x2 = mesh_net.meshes[i].p_tar.p.x();
			double y2 = mesh_net.meshes[i].p_tar.p.y();
			if (x1<min_x){min_x=x1;}
			if (x1>max_x){max_x=x1;}
			if (y1<min_y){min_y=y1;}
			if (y1>max_y){max_y=y1;}
			if (x2<min_x){min_x=x2;}
			if (x2>max_x){max_x=x2;}
			if (y2<min_y){min_y=y2;}
			if (y2>max_y){max_y=y2;}
		}
		cout << "DFN extents: min(" << min_x << "," << min_y << "), max("
		     << max_x << "," << max_y << ")" << endl;
		cout << "DFN vs domain x-gap: left " << (min_x - domain.min_pt.i)
		     << ", right " << (domain.max_pt.i - max_x) << endl;
	};
	pair<int,int> borders_before = count_border_meshes(init_net_mesh);
	cout << "DFN meshes before backbone: " << mesh_count_before_backbone
	     << ", left-border meshes: " << borders_before.first
	     << ", right-border meshes: " << borders_before.second
	     << ", meshes inside domain: " << count_in_domain_meshes(init_net_mesh) << endl;
	report_dfn_extents(init_net_mesh);
	// Preserve the generated network before backbone filtering for inspection.
	{
		ofstream raw_dfn_output((param.code_path+"/Output/DFN_raw.txt").c_str(),ofstream::out);
		for (size_t i=0;i<init_net_mesh.meshes.size();i++){
			raw_dfn_output << init_net_mesh.meshes[i].p_ori.p.x() << " "
			               << init_net_mesh.meshes[i].p_ori.p.y() << " "
			               << init_net_mesh.meshes[i].p_tar.p.x() << " "
			               << init_net_mesh.meshes[i].p_tar.p.y() << " "
			               << init_net_mesh.meshes[i].aperture << endl;
		}
		raw_dfn_output.close();
	}
	if (mesh_count_before_backbone==0){
		cout << "WARNING in PERFORM.cpp: DFN mesh count is 0 after file read; check file format and coordinates" << endl;
	}
	
	if (backbone){
	    cout << "Fluid flow computation in the fracture network" << endl;
		double velocity_threshold=EPSILON;
		net_mesh=net_mesh.return_backbone(param,velocity_threshold);
		init_net_mesh=net_mesh;
		if (net_mesh.meshes.size()>0){
			FractureMeshNumbering(net_mesh);
			UpdateCptFract(net_mesh); // update the map of fractures on the backbone
		}
		else if (mesh_count_before_backbone>0){
			cout << "WARNING in PERFORM.cpp: backbone filtering removed all meshes; DFN connectivity lost" << endl;
		}
		pair<int,int> borders_after = count_border_meshes(net_mesh);
		cout << "DFN meshes after backbone: " << net_mesh.meshes.size()
		     << ", left-border meshes: " << borders_after.first
		     << ", right-border meshes: " << borders_after.second
		     << ", meshes inside domain: " << count_in_domain_meshes(net_mesh) << endl;
		report_dfn_extents(net_mesh);
	}
	auto connectivity_end = std::chrono::steady_clock::now();
	cout << "Timing: initial connectivity+flow = "
	     << std::chrono::duration<double>(connectivity_end-connectivity_start).count() << " s" << endl;

	// 3. Transport modeling
	cout << "2. Transport process computation - Particle displacement" << endl;
	RngStream_a rng_tracker;
	Transport transport(rng_tracker,net_mesh,domain,param);
	if (option!=0&&net_mesh.meshes.size()>0){
		map<int,double> arrival_times;
		if (transport.Particles_Transport(arrival_times,option_injection)){
		        cout << "3. Result post-processing and writting" << endl;
			int arrived = arrival_times.size();
			int total = param.nb_part;
			double min_arrival = 0.0;
			double max_arrival = 0.0;
			double mean_arrival = 0.0;
			if (!arrival_times.empty()){
				min_arrival = arrival_times.begin()->second;
				max_arrival = arrival_times.begin()->second;
				for (map<int,double>::iterator it=arrival_times.begin(); it!=arrival_times.end(); it++){
					double t = it->second;
					if (t<min_arrival){min_arrival = t;}
					if (t>max_arrival){max_arrival = t;}
					mean_arrival += t;
				}
				mean_arrival /= (double)arrival_times.size();
			}
			cout << "Diagnostics: arrived=" << arrived
			     << " total=" << total
			     << " not_arrived=" << (total-arrived)
			     << " mean_t=" << mean_arrival
			     << " min_t=" << min_arrival
			     << " max_t=" << max_arrival << endl;
			int sample = 0;
			for (map<int,double>::iterator it=arrival_times.begin(); it!=arrival_times.end() && sample<5; it++,sample++){
				cout << "Diagnostics sample arrival: particle=" << it->first
				     << " t=" << it->second << endl;
			}
			transport.PrintTravelTimeSamplesAtOutputIntervals(5);
			// 5. Results post-processing
			auto results_start = std::chrono::steady_clock::now();
			Results results(arrival_times,param);
			results.post_processing();
			results.writing(param.code_path);
			transport.WritePositionSnapshotsCSV(param.code_path);
			auto results_end = std::chrono::steady_clock::now();
			cout << "Timing: results output = "
			     << std::chrono::duration<double>(results_end-results_start).count() << " s" << endl;
		}
	}
	
        // 4. Display fracture network
	if (net_mesh.meshes.size()==0){cout << "; WARNING in PERFORM.cpp: not connected backbone";}
	else{
		//net_mesh.print_DFN_in_file(param);
		transport.net_mesh_modified.print_DFN_in_file(param);  // modified by DR on 2025/12/11
		transport.net_mesh_modified.print_DFN_in_file_with_aperture_delta(param,init_net_mesh);
		init_net_mesh.print_DFN_in_file_initial(param);
		//DisplayResults(argc,argv,net_mesh);
	}
	cout << "CPU Time = " << (clock()-t_begin_total)/CLOCKS_PER_SEC << endl;

	return 0;
}
