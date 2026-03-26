/*
 * Transport.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "Transport.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <set>
#include "../Domain_Definition/DFNComputation.h"
using namespace std;

static const double VELOCITY_APERTURE_ZERO_THRESHOLD = 1e-10;

// Modified by Wenyu on 2026/1/12; check if a mesh has connected meshes at a given node
static bool HasConnectedMeshAtNode(const NetworkMeshes& net_mesh,int mesh_index,int node_index){
	FractureMesh mesh = net_mesh.return_mesh(mesh_index);
	for (size_t i=0;i<mesh.neigh_meshes.size();i++){
		FractureMesh neigh = net_mesh.return_mesh(mesh.neigh_meshes[i]);
		if (neigh.p_ori.index==node_index || neigh.p_tar.index==node_index){
			return true;
		}
	}
	return false;
}
// Modified by Wenyu on 2026/1/12; rebuild the network that removing meshes with zero aperture 
// (Directly modified by variable not from the functions)
static void RebuildNetworkRemovingZeroAperture(NetworkMeshes& net_mesh,
                                               std::map<int,int>& old_to_new,
                                               std::map<int,int>& node_old_to_new,
                                               std::vector<int>& removed_meshes,
                                               double& rebuild_seconds){
	auto t_start = std::chrono::steady_clock::now();
	std::vector<FractureMesh> new_meshes;
	std::vector<int> old_mesh_indices;
	new_meshes.reserve(net_mesh.meshes.size());
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		if (net_mesh.meshes[i].aperture>0.0){ // keep the mesh above zero 
			old_mesh_indices.push_back(net_mesh.meshes[i].mesh_index);
			new_meshes.push_back(net_mesh.meshes[i]);
		}
		else{
			removed_meshes.push_back(net_mesh.meshes[i].mesh_index);
		}
	}
	net_mesh.meshes=new_meshes;
	// Rebuild node indices and maps; Modified by Wenyu on 2026/1/13
	std::set<int> node_set;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		node_set.insert(net_mesh.meshes[i].p_ori.index);
		node_set.insert(net_mesh.meshes[i].p_tar.index);
	}
	std::map<int,int> node_map;
	int node_idx = 0;
	for (std::set<int>::iterator it=node_set.begin(); it!=node_set.end(); it++){
		node_map[*it]=node_idx++;
	}
	node_old_to_new = node_map;
	net_mesh.nodes_map.clear();
	net_mesh.border_map.clear();
	net_mesh.inter_list.clear();
	std::map<int,int> node_degree;
	for (size_t i=0;i<net_mesh.meshes.size();i++){ // update node indices in meshes
		FractureMesh& mesh = net_mesh.meshes[i];
		mesh.p_ori.index = node_map[mesh.p_ori.index];
		mesh.p_tar.index = node_map[mesh.p_tar.index];
		net_mesh.nodes_map[mesh.p_ori.index]=mesh.p_ori.p;
		net_mesh.nodes_map[mesh.p_tar.index]=mesh.p_tar.p;
		node_degree[mesh.p_ori.index] += 1;
		node_degree[mesh.p_tar.index] += 1;
	}
	// Rebuild border map
	for (std::map<int,CgalPoint2D>::iterator it=net_mesh.nodes_map.begin(); it!=net_mesh.nodes_map.end(); it++){
		std::string border=net_mesh.domain.ReturnBorder(it->second);
		if (border!=NO_BORDER){net_mesh.border_map[it->first]=border;}
	}
	// Rebuild inter_list
	for (std::map<int,int>::iterator it=node_degree.begin(); it!=node_degree.end(); it++){
		if (it->second>1){net_mesh.inter_list.insert(it->first);}
	}
	// Rebuild mesh indices and neigh_meshes
	net_mesh.cpt_inter = node_idx-1;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		net_mesh.meshes[i].neigh_meshes.clear();
	}
	FractureMeshNumbering(net_mesh);
	UpdateCptFract(net_mesh); // to update cpt_fract
	// Rebuild Nb_passed and last_t_update
	for (size_t i=0;i<net_mesh.meshes.size() && i<old_mesh_indices.size();i++){
		old_to_new[old_mesh_indices[i]] = (int)net_mesh.meshes[i].mesh_index;
	}
	if (!net_mesh.Nb_passed.empty()){
		std::vector<int> new_nb(net_mesh.meshes.size(),0);
		for (std::map<int,int>::iterator it=old_to_new.begin(); it!=old_to_new.end(); it++){
			if (it->first>=0 && it->first<(int)net_mesh.Nb_passed.size()){
				new_nb[it->second]=net_mesh.Nb_passed[it->first];
			}
		}
		net_mesh.Nb_passed=new_nb;
	}
	// Rebuild last_t_update
	if (!net_mesh.last_t_update.empty()){
		std::vector<double> new_last(net_mesh.meshes.size(),0.0);
		for (std::map<int,int>::iterator it=old_to_new.begin(); it!=old_to_new.end(); it++){
			if (it->first>=0 && it->first<(int)net_mesh.last_t_update.size()){
				new_last[it->second]=net_mesh.last_t_update[it->first];
			}
		}
		net_mesh.last_t_update=new_last;
	}
	auto t_end = std::chrono::steady_clock::now();
	rebuild_seconds = std::chrono::duration<double>(t_end-t_start).count();
}
static void RebuildNodeNumbering(NetworkMeshes& net_mesh,std::map<int,int>& old_to_new_mesh){
	std::vector<int> old_mesh_indices;
	old_mesh_indices.reserve(net_mesh.meshes.size());
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		old_mesh_indices.push_back(net_mesh.meshes[i].mesh_index);
	}
	std::set<int> node_set;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		node_set.insert(net_mesh.meshes[i].p_ori.index);
		node_set.insert(net_mesh.meshes[i].p_tar.index);
	}
	std::map<int,int> node_map;
	int node_idx = 0;
	for (std::set<int>::iterator it=node_set.begin(); it!=node_set.end(); it++){
		node_map[*it]=node_idx++;
	}
	net_mesh.nodes_map.clear();
	net_mesh.border_map.clear();
	net_mesh.inter_list.clear();
	std::map<int,int> node_degree;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		FractureMesh& mesh = net_mesh.meshes[i];
		mesh.p_ori.index = node_map[mesh.p_ori.index];
		mesh.p_tar.index = node_map[mesh.p_tar.index];
		net_mesh.nodes_map[mesh.p_ori.index]=mesh.p_ori.p;
		net_mesh.nodes_map[mesh.p_tar.index]=mesh.p_tar.p;
		node_degree[mesh.p_ori.index] += 1;
		node_degree[mesh.p_tar.index] += 1;
	}
	for (std::map<int,CgalPoint2D>::iterator it=net_mesh.nodes_map.begin(); it!=net_mesh.nodes_map.end(); it++){
		std::string border=net_mesh.domain.ReturnBorder(it->second);
		if (border!=NO_BORDER){net_mesh.border_map[it->first]=border;}
	}
	for (std::map<int,int>::iterator it=node_degree.begin(); it!=node_degree.end(); it++){
		if (it->second>1){net_mesh.inter_list.insert(it->first);}
	}
	net_mesh.cpt_inter = node_idx-1;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		net_mesh.meshes[i].neigh_meshes.clear();
	}
	FractureMeshNumbering(net_mesh);
	UpdateCptFract(net_mesh);
	old_to_new_mesh.clear();
	for (size_t i=0;i<net_mesh.meshes.size() && i<old_mesh_indices.size();i++){
		old_to_new_mesh[old_mesh_indices[i]] = net_mesh.meshes[i].mesh_index;
	}
	if (!net_mesh.Nb_passed.empty()){
		std::vector<int> new_nb(net_mesh.meshes.size(),0);
		for (std::map<int,int>::iterator it=old_to_new_mesh.begin(); it!=old_to_new_mesh.end(); it++){
			if (it->first>=0 && it->first<(int)net_mesh.Nb_passed.size()){
				new_nb[it->second]=net_mesh.Nb_passed[it->first];
			}
		}
		net_mesh.Nb_passed=new_nb;
	}
	if (!net_mesh.last_t_update.empty()){
		std::vector<double> new_last(net_mesh.meshes.size(),0.0);
		for (std::map<int,int>::iterator it=old_to_new_mesh.begin(); it!=old_to_new_mesh.end(); it++){
			if (it->first>=0 && it->first<(int)net_mesh.last_t_update.size()){
				new_last[it->second]=net_mesh.last_t_update[it->first];
			}
		}
		net_mesh.last_t_update=new_last;
	}
}
// Modified by Wenyu on 2026/1/13; return the closest endpoint of a mesh to a reference point (NOT USED)
static pointcpp<double> ClosestMeshEndpoint(const FractureMesh& mesh,const pointcpp<double>& ref){
	pointcpp<double> p1(mesh.p_ori.p);
	pointcpp<double> p2(mesh.p_tar.p);
	double d1 = ref.Point_Distance(p1);
	double d2 = ref.Point_Distance(p2);
	return (d1<=d2)? p1 : p2;
}
static bool PointsClose(const pointcpp<double>& a,const pointcpp<double>& b){
	return a.Point_Distance(b)<=EPSILON;
}
static bool SameSegment(const FractureMesh& a,const FractureMesh& b){
	pointcpp<double> a1(a.p_ori.p), a2(a.p_tar.p);
	pointcpp<double> b1(b.p_ori.p), b2(b.p_tar.p);
	if (PointsClose(a1,b1) && PointsClose(a2,b2)){return true;}
	if (PointsClose(a1,b2) && PointsClose(a2,b1)){return true;}
	return false;
}
static int FindMatchingMeshIndex(const NetworkMeshes& new_mesh,const FractureMesh& old_mesh){
	for (size_t i=0;i<new_mesh.meshes.size();i++){
		const FractureMesh& candidate = new_mesh.meshes[i];
		if (candidate.fracture_index!=old_mesh.fracture_index){continue;}
		if (SameSegment(candidate,old_mesh)){return candidate.mesh_index;}
	}
	return -1;
}
static void ResetParticleForInlet(Particle& pa,double t_reset){
	pa.t = t_reset;
	pa.t_injection = t_reset;
	pa.mesh_index = -1;
	pa.M = pointcpp<double>(NOT_DEFINED,NOT_DEFINED);
	pa.mesh_history.clear();
	pa.prev_mesh_index = -1;
	pa.intersection_history.clear();
	pa.L_in_fract = 0.0;
	pa.t_in_fract = 0.0;
	pa.t_in_fract_prev = 0.0;
}
static void WriteDFNSnapshot(const NetworkMeshes& net_mesh,const std::string& output_path,
                             int file_index,double time){
	std::ostringstream name_stream;
	name_stream << output_path << "/Output/DFN_step"
	            << std::setw(6) << std::setfill('0') << file_index << ".txt";
	std::string file_name = name_stream.str();
	std::ofstream output(file_name.c_str(), std::ofstream::out);
	if (!output.is_open()){
		cout << "WARNING in WriteDFNSnapshot (Transport.cpp): cannot open file " << file_name << endl;
		return;
	}
	output << "# time=" << time << endl;
	output << "x1 y1 x2 y2 aperture velocity" << endl;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		const FractureMesh& mesh = net_mesh.meshes[i];
		output << mesh.p_ori.p.x() << " " << mesh.p_ori.p.y() << " "
		       << mesh.p_tar.p.x() << " " << mesh.p_tar.p.y() << " "
		       << mesh.aperture << " " << mesh.velocity << endl;
	}
	output.close();
}
// Modified by Wenyu on 2026/1/19; remove particle from the system (USED)
static void RemoveParticleFromSystem(Particle& pa,double t_stamp,const std::string& reason){
	cout << "Particle removed: particle=" << pa.no
	     << " mesh_index=" << pa.mesh_index
	     << " t=" << t_stamp
	     << " reason=" << reason << endl;
	pa.t = -1;
	pa.mesh_index = -1;
	pa.M = pointcpp<double>(NOT_DEFINED,NOT_DEFINED);
	pa.mesh_history.clear();
	pa.prev_mesh_index = -1;
	pa.intersection_history.clear();
	pa.L_in_fract = 0.0;
	pa.t_in_fract = 0.0;
	pa.t_in_fract_prev = 0.0;
}
// Modified by Wenyu on 2026/1/15; reassign particle to the nearest connected mesh endpoint (USED)
static bool ReassignToNearestConnectedEndpoint(Particle& pa,const NetworkMeshes& net_mesh,
                                               Domain& domain,
                                               const std::set<int>& connected_nodes,
                                               double t_target,bool update_time){
	int old_mesh = pa.mesh_index;
	// check if position is valid， reassin to inlet if not (NaN)
	if (!std::isfinite(pa.M.i) || !std::isfinite(pa.M.j)){
		cout << "WARNING: particle " << pa.no
		     << " has invalid position, reassign to inlet from mesh " << old_mesh
		     << " t=" << pa.t << endl;
		ResetParticleForInlet(pa,update_time ? t_target : pa.t);
		return true;
	}
	double best_dist = 1e100;
	int best_mesh = -1;
	pointcpp<double> best_pt;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		const FractureMesh& mesh = net_mesh.meshes[i];
		if (mesh.mesh_index==old_mesh){
			continue;
		}
		if (connected_nodes.find(mesh.p_ori.index)==connected_nodes.end() ||
		    connected_nodes.find(mesh.p_tar.index)==connected_nodes.end()){
			continue;
		}
		pointcpp<double> p1(mesh.p_ori.p);
		pointcpp<double> p2(mesh.p_tar.p);
		bool p1_ok = std::isfinite(p1.i) && std::isfinite(p1.j);
		bool p2_ok = std::isfinite(p2.i) && std::isfinite(p2.j);
		if (p1_ok){
			double d1 = pa.M.Point_Distance(p1);
			if (d1<best_dist){
				best_dist = d1;
				best_mesh = mesh.mesh_index;
				best_pt = p1;
			}
		}
		if (p2_ok){
			double d2 = pa.M.Point_Distance(p2);
			if (d2<best_dist){
				best_dist = d2;
				best_mesh = mesh.mesh_index;
				best_pt = p2;
			}
		}
	}
	if (best_mesh<0){
		cout << "WARNING: particle " << pa.no
		     << " removed (no connected mesh) from mesh " << old_mesh
		     << " t=" << pa.t << endl;
		pa.t = -1;
		pa.mesh_history.clear();
		pa.prev_mesh_index = -1;
		pa.intersection_history.clear();
		return false;
	}
	if (domain.on_input_limit(best_pt)){
		ResetParticleForInlet(pa,update_time ? t_target : pa.t);
		return true;
	}
	pa.mesh_index = best_mesh;
	pa.M = best_pt;
	pa.mesh_history.clear();
	pa.prev_mesh_index = -1;
	pa.intersection_history.clear();
	if (update_time){
		double dt_remain = t_target - pa.t;
		if (dt_remain>0.0){
			pa.t = t_target;
			pa.t_in_fract += dt_remain;
		}
		cout << "Particle reassigned to nearest connected mesh endpoint: particle=" << pa.no
		     << " from mesh " << old_mesh
		     << " to mesh " << best_mesh
		     << " step_t=" << t_target << endl;
	}
	else{
		cout << "Particle reassigned to nearest connected mesh endpoint: particle=" << pa.no
		     << " from mesh " << old_mesh
		     << " to mesh " << best_mesh
		     << " t=" << pa.t << endl;
	}
	return true;
}
// Modified by Wenyu on 2026/1/12; Normalize mesh directions so that all velocities are positive
static void NormalizeMeshDirections(NetworkMeshes& net_mesh){
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		if (net_mesh.meshes[i].velocity<0.0){
			FractureMesh tmp = net_mesh.meshes[i];
			net_mesh.meshes[i].p_ori = tmp.p_tar;
			net_mesh.meshes[i].p_tar = tmp.p_ori;
			net_mesh.meshes[i].velocity = -tmp.velocity;
		}
	}
}
// Modified by Wenyu on 2026/1/13; set aperture to zero if velocity is below threshold
static bool ZeroApertureIfLowVelocity(NetworkMeshes& net_mesh){
	bool changed = false;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		if (fabs(net_mesh.meshes[i].velocity) <= VELOCITY_APERTURE_ZERO_THRESHOLD &&
		    net_mesh.meshes[i].aperture>0.0){
			net_mesh.meshes[i].aperture = 0.0;
			changed = true;
		}
	}
	return changed;
}
// Modified by Wenyu on 2026/1/14; collect initial positions for injection
static void CollectInputPositions(const NetworkMeshes& net_mesh,Domain& domain,
                                  std::map<int,pointcpp<double> >& init_pos){
	init_pos.clear();
	pointcpp<double> pt1, pt2;
	for (size_t i=0;i<net_mesh.meshes.size();i++){
		const FractureMesh& fract_mesh = net_mesh.meshes[i];
		if (fract_mesh.aperture<=0.0){continue;}
		pt1 = pointcpp<double>(fract_mesh.p_ori.p);
		pt2 = pointcpp<double>(fract_mesh.p_tar.p);
		if (domain.on_input_limit(pt1)){
			init_pos[fract_mesh.mesh_index] = pt1;
		}
		else if (domain.on_input_limit(pt2)){
			init_pos[fract_mesh.mesh_index] = pt2;
		}
	}
}
// Modified by Wenyu on 2026/1/14; compute injection counts based on weights
static bool ComputeInjectionCounts(int total,const std::vector<double>& weights,
                                   std::vector<int>& counts){
	int nb_pos = weights.size();
	counts.assign(nb_pos,0);
	if (total<=0 || nb_pos==0){return false;}
	double sum_w = 0.0;
	for (int i=0;i<nb_pos;i++){sum_w += weights[i];}
	if (sum_w<=0.0){return false;}
	std::vector<double> frac(nb_pos,0.0);
	int base_sum = 0;
	for (int i=0;i<nb_pos;i++){
		double exact = total*weights[i]/sum_w;
		int base = (int)floor(exact);
		counts[i]=base;
		base_sum+=base;
		frac[i]=exact-base;
	}
	int leftover = total-base_sum;
	std::vector<int> order(nb_pos,0);
	for (int i=0;i<nb_pos;i++){order[i]=i;}
	sort(order.begin(),order.end(),[&](int a,int b){
		if (frac[a]==frac[b]){return a<b;}
		return frac[a]>frac[b];
	});
	for (int k=0;k<leftover;k++){
		counts[order[k%nb_pos]]++;
	}
	return true;
}
// Modified by Wenyu on 2026/1/14; print injection counts
static void PrintInjectionCounts(const std::string& label,
                                 const std::vector<int>& mesh_ids,
                                 const std::vector<int>& counts){
	cout << label << " (mesh_index -> count):" << endl;
	for (size_t i=0;i<mesh_ids.size();i++){
		cout << mesh_ids[i] << " -> " << counts[i] << endl;
	}
}
// Modified by Wenyu on 2026/1/14; inject particles at initial positions based on weights
Transport::Transport(RngStream_a rng_tracker_,NetworkMeshes net_mesh_,Domain domain_,Parameters param){
	rng_tracker = rng_tracker_;
	net_mesh = net_mesh_;
	domain = domain_;
	param_full = param;
	phys_param = PhysicsParam(param.Dm,param.porosity);
	num_param = NumericalParam(param.nb_part,param.proba_transfer,param.simu_option,param.t_max,param.t_injection,param.output_interval,param.reaction_dt);
	net_mesh_modified = net_mesh; // added by DR on 2025/12/11
	net_mesh_modified.Nb_passed.assign(net_mesh_modified.meshes.size(), 0);// added by Wenyu on 2025/12/30
	net_mesh_modified.last_t_update.assign(net_mesh_modified.meshes.size(), 0.0);// added by Wenyu on 2025/12/30

};

// injection of particle on the left border of the domain
vector<Particle> Transport::Particles_Injection(int option_injection){
	vector<Particle> part_vect(this->num_param.nb_part);

	int nb_part_tot = part_vect.size();
	double t_injection = num_param.t_injection;
	double dt_inject = 0.0;
	if (t_injection>0.0 && nb_part_tot>1){
		dt_inject = t_injection/(double)(nb_part_tot-1);
	}
	if (t_injection>0.0){
		cout << "Continuous injection enabled: t_injection=" << t_injection
		     << ", dt_inject=" << dt_inject << endl;
	}

	for (int i=0;i<nb_part_tot;i++){
		part_vect[i].M = pointcpp<double>(NOT_DEFINED,NOT_DEFINED);
		part_vect[i].mesh_index = -1;
		part_vect[i].prev_mesh_index = -1;
		part_vect[i].mesh_history.clear();
		part_vect[i].t = (t_injection>0.0)? (i*dt_inject) : 0.0;
		part_vect[i].t_injection = part_vect[i].t;
		part_vect[i].no = i;
		part_vect[i].L_in_fract = 0.0;
		part_vect[i].t_in_fract = 0.0;
		part_vect[i].t_in_fract_prev = 0.0;
	}
	return part_vect;
}

// particle displacement considering the surrounding fractures
bool Transport::finite_matrix_displacement(Particle & pa){
	// 0. Variables
	// 0.1. Inputs
	FractureMesh current_mesh = net_mesh.return_mesh(pa.mesh_index);
	pointcpp<double> M_init = pa.M, M_init_proj, M_out=pointcpp<double>(NOT_DEFINED,NOT_DEFINED), M_proj,M_init_in;pair<pointcpp<double>,pointcpp<double> > seg_current;
	double porosity = phys_param.porosity, Dm = phys_param.Dm, transfer_time,particle_time,current_abs_velocity;
	double t_max_advec,t_advec,t_diff;
	int transfer_mesh=-1;

	Projection_Map projections;
	int nb_iter=0,nb_max_iter=1e4;

	// 1. Particle displacement until an intersection
	while (!M_init.identic(M_out, EPSILON)&&domain.IsInDomain(M_init.CgalPoint())&&pa.t<num_param.t_max&&nb_iter<nb_max_iter){	// MODIF DR 2020-11-03
		nb_iter++;
		// 1.1. Determination of the characteristic times (depending on fracture discretization)
		// 1.1.1. Determination of the fracture extremity (in the sense of the flow)
		current_mesh.distance_from_M_to_extremity(M_init, M_out);
		t_max_advec = current_mesh.Mesh_Time_Advection(M_init,M_out);
		// 1.1.2. Determination of the neighboring fractures
		M_init_proj = current_mesh.Mesh_Scale_Advection(M_init,0.5*t_max_advec);	//small displacement to get the "good" meshes of projection
		
		//Projection_Map projections = Orthogonal_Projection_M(M_init_proj.CgalPoint(), M_out.CgalPoint());
		if (!exist_in_map(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections)){
			/*cout << "Projection 1" << endl;
                        cout << "M_init = ";print(M_init.CgalPoint());
			cout << "M_init_proj = ";print(M_init_proj.CgalPoint());
                        cout << "M_out = ";print(M_out.CgalPoint());*/
			//projections=Orthogonal_Projection_M(M_init_proj.CgalPoint(), M_out.CgalPoint());
			Projection_Map projections_new;
			if (!Orthogonal_Projection_M(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections_new)){return false;}
			projections=projections_new;
			proj_map_global[make_pair(M_init_proj.CgalPoint(), M_out.CgalPoint())]=projections;
		}
		// 1.1.3. Determination of the corresponding advection and diffusion times (depending on projected distance and transfer_proba_max)
		t_advec = current_mesh.Advection_Time_Computation(M_init,projections,porosity,Dm,num_param.proba_transfer,t_max_advec);
		if (!current_mesh.Get_Total_Time_From_Advec_Time1D(t_diff,t_advec,Dm,porosity,rng_tracker,pa.t_in_fract)){return false;};
		t_diff=t_diff-t_advec;

		//t_diff = current_mesh.Get_Total_Time_From_Advec_Time(t_advec,Dm,porosity,rng_tracker,previous_time,Initial_Position_Particles[pa.no],M_out)-t_advec;
		// 1.2. Corresponding displacements
		// 1.2.1. Displacement by advection (and update of the orthogonal projections)
		M_init_in=M_init;
		M_init = current_mesh.Mesh_Scale_Advection(M_init, t_advec);
		seg_current=make_pair(M_init_in,M_init);current_abs_velocity=abs(current_mesh.velocity);
		if (M_init.identic(M_out, EPSILON)){M_init_proj = current_mesh.Mesh_Scale_Advection(M_init, t_advec/10, true);}
		else{M_init_proj = M_init;}
		//projections = Orthogonal_Projection_M(M_init_proj.CgalPoint(), M_out.CgalPoint());
		if (!exist_in_map(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections)){
			//projections=Orthogonal_Projection_M(M_init_proj.CgalPoint(),M_out.CgalPoint());
			Projection_Map projections_new;
			if (!Orthogonal_Projection_M(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections_new)){return false;}
			projections=projections_new;
			proj_map_global[make_pair(M_init_proj.CgalPoint(), M_out.CgalPoint())]=projections;
		}

		// 1.2.2. Displacement by diffusion
		// if the particle transfers to a neighbor fracture by diffusion
		if (projections.size()>0 && Transfer_Probability_and_Transfer_Time(projections,t_diff,transfer_mesh,transfer_time,M_proj)){
			particle_time = t_advec+transfer_time;	// particle time updating
			// standard case: the particle transfers to a new fracture
			if (transfer_mesh!=-1){
				current_mesh = net_mesh.return_mesh(transfer_mesh);	// mesh current updating
				pa.mesh_index = transfer_mesh;
				M_init=M_proj;// particle position updating
				// distance from the new current position of the particle to the origin of the fracture in which the particle just entered by diffusion
				pa.L_in_fract=M_init.Point_Distance(net_mesh.fractures[current_mesh.fracture_index].first.p);
			}
			// boundary case: the particle reflects on a domain border
			else{cout << "WARNING in finite_matrix_displacement (Transport.cpp): transfer mesh not defined" << endl;}
		}
		// if the particle does not transfer (go back to the initial fracture after diffusion in the matrix)
		else{// update particle time and covered distance in the current fracture
			particle_time = t_advec+t_diff;
			pa.t_in_fract+=particle_time;
			pa.L_in_fract+=fabs(current_mesh.velocity)*t_advec;
		}
		// 1.2.3. Updating of particle and positions characteristics
		pa.t += particle_time;pa.M = M_init;
		current_mesh.distance_from_M_to_extremity(M_init, M_out);
		Segment_Particle_Time[seg_current].insert(pa.t); // when the particle reaches another position
		Segment_Velocities[seg_current]=current_abs_velocity;
	}

	if (nb_iter==nb_max_iter||pa.t>=num_param.t_max){cout << "Particle stopped because of nb_max_iter or t_max, simulation stopped"; pa.t=-1; return false;}

	return true;
}
// Modified by Wenyu on 2026/1/8: Finite matrix displacement with distance tracking (to line 352)
bool Transport::finite_matrix_displacement_step(Particle & pa, double t_target, std::map<int,double>& moved_distances){
	// 0. Variables
	FractureMesh current_mesh = net_mesh.return_mesh(pa.mesh_index);
	static std::set<int> warned_zero_velocity;
	if ((current_mesh.velocity==0.0 || current_mesh.aperture==0.0) &&
	    warned_zero_velocity.insert(pa.mesh_index).second){
		cout << "Attention: in finite_matrix_displacement_step aperture closed == 0: mesh_index=" << pa.mesh_index
		     << " velocity=" << current_mesh.velocity
		     << " aperture=" << current_mesh.aperture << endl;
	}
	// If velocity is zero, just advance time to t_target; not stop simulation, update the flow field
	if (fabs(current_mesh.velocity==0.0)){
		pointcpp<double> M_end(current_mesh.p_tar.p);
		pointcpp<double> M_start(current_mesh.p_ori.p);
		if (!HasConnectedMeshAtNode(net_mesh,pa.mesh_index,current_mesh.p_tar.index)){
			pa.M = M_start;
		}
		else{
			pa.M = M_end;
		}
		double dt_remain = t_target - pa.t;
		if (dt_remain>0.0){pa.t_in_fract+=dt_remain;}
		pa.t = t_target;
		return true;
	}
	pointcpp<double> M_init = pa.M, M_init_proj, M_out=pointcpp<double>(NOT_DEFINED,NOT_DEFINED), M_proj,M_init_in;pair<pointcpp<double>,pointcpp<double> > seg_current;
	double porosity = phys_param.porosity, Dm = phys_param.Dm, transfer_time,particle_time,current_abs_velocity;
	double t_max_advec,t_advec,t_diff;
	int transfer_mesh=-1;

	Projection_Map projections;
	int nb_iter=0,nb_max_iter=1e4;

	while (!M_init.identic(M_out, EPSILON)&&domain.IsInDomain(M_init.CgalPoint())&&pa.t<t_target&&nb_iter<nb_max_iter){
		nb_iter++;
		current_mesh.distance_from_M_to_extremity(M_init, M_out);
		t_max_advec = current_mesh.Mesh_Time_Advection(M_init,M_out);
		M_init_proj = current_mesh.Mesh_Scale_Advection(M_init,0.5*t_max_advec);

		if (!exist_in_map(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections)){
			Projection_Map projections_new;
			if (!Orthogonal_Projection_M(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections_new)){return false;}
			projections=projections_new;
			proj_map_global[make_pair(M_init_proj.CgalPoint(), M_out.CgalPoint())]=projections;
		}
		t_advec = current_mesh.Advection_Time_Computation(M_init,projections,porosity,Dm,num_param.proba_transfer,t_max_advec);
		if (!current_mesh.Get_Total_Time_From_Advec_Time1D(t_diff,t_advec,Dm,porosity,rng_tracker,pa.t_in_fract)){return false;}
		t_diff=t_diff-t_advec;

		M_init_in=M_init;
		M_init = current_mesh.Mesh_Scale_Advection(M_init, t_advec);
		seg_current=make_pair(M_init_in,M_init);current_abs_velocity=abs(current_mesh.velocity);
		if (M_init.identic(M_out, EPSILON)){M_init_proj = current_mesh.Mesh_Scale_Advection(M_init, t_advec/10, true);}
		else{M_init_proj = M_init;}
		if (!exist_in_map(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections)){
			Projection_Map projections_new;
			if (!Orthogonal_Projection_M(M_init_proj.CgalPoint(),M_out.CgalPoint(),projections_new)){return false;}
			projections=projections_new;
			proj_map_global[make_pair(M_init_proj.CgalPoint(), M_out.CgalPoint())]=projections;
		}

		bool transferred = false;
		if (projections.size()>0 && Transfer_Probability_and_Transfer_Time(projections,t_diff,transfer_mesh,transfer_time,M_proj)){
			particle_time = t_advec+transfer_time;
			transferred = true;
		}
		else{
			particle_time = t_advec+t_diff;
		}

		if (pa.t + particle_time > t_target){
			double dt_remain = t_target - pa.t;
			double adv_time_used = std::min(dt_remain, t_advec);
			M_init = current_mesh.Mesh_Scale_Advection(M_init_in, adv_time_used);
			pa.t = t_target;
			pa.M = M_init;
			pa.t_in_fract+=dt_remain;
			pa.L_in_fract+=fabs(current_mesh.velocity)*adv_time_used;
			moved_distances[pa.mesh_index]+=fabs(current_mesh.velocity)*adv_time_used;
			return true;
		}

		moved_distances[pa.mesh_index]+=fabs(current_mesh.velocity)*t_advec;
		if (transferred && transfer_mesh!=-1){
			current_mesh = net_mesh.return_mesh(transfer_mesh);
			pa.mesh_index = transfer_mesh;
			M_init=M_proj;
			pa.L_in_fract=M_init.Point_Distance(net_mesh.fractures[current_mesh.fracture_index].first.p);
		}
		else if (transferred && transfer_mesh==-1){
			cout << "WARNING in finite_matrix_displacement_step (Transport.cpp): transfer mesh not defined" << endl;
		}
		else{
			pa.t_in_fract+=particle_time;
			pa.L_in_fract+=fabs(current_mesh.velocity)*t_advec;
		}

		pa.t += particle_time;pa.M = M_init;
		current_mesh.distance_from_M_to_extremity(M_init, M_out);
		Segment_Particle_Time[seg_current].insert(pa.t);
		Segment_Velocities[seg_current]=current_abs_velocity;
	}

	if (nb_iter==nb_max_iter||pa.t>=num_param.t_max){cout << "Particle stopped because of nb_max_iter or t_max, simulation stopped"; pa.t=-1; return false;}

	return true;
}

bool Transport::exist_in_map(pointcpp<double> pt1,pointcpp<double> pt2, Projection_Map & projections){
	Projection_Map_Global::iterator it=proj_map_global.find(make_pair(pt1,pt2));
	if (it!=proj_map_global.end()){
		projections=it->second;
		return true;
	}
	return false;
}

// particle displacement considering an infinite surrounding matrix
bool Transport::infinite_matrix_displacement(Particle & pa){
	// 0. Variables
	FractureMesh current_mesh = net_mesh.return_mesh(pa.mesh_index);

	if (current_mesh.velocity<0&&pa.M.identic(pointcpp<double>(current_mesh.p_ori.p),EPSILON)){
		cout << "WARNING in infinite_matrix_displacement: flow velocity negative from origin to target" << endl;exit(0);
	}

	pointcpp<double> M_init = pa.M, M_out;
	// 1. Displacement from the first to the second extremity of the mesh
	// definition of M_out (extremity opposite to M_init)
	current_mesh.distance_from_M_to_extremity(M_init, M_out);
	double advection_time = current_mesh.Mesh_Time_Advection(M_init, M_out);
	double particle_time;
	if (pa.t==-1){pa.M = M_out;}
	else{
		//if (!current_mesh.Get_Total_Time_From_Advec_Time(particle_time,advection_time,phys_param.Dm,phys_param.porosity,rng_tracker,pa.t_in_fract,pa.L_in_fract,num_param.Nt_arr,Delta_CDF_Map,pa.no)){return false;}
		if (!current_mesh.Get_Total_Time_From_Advec_Time1D(particle_time,advection_time,phys_param.Dm,phys_param.porosity,rng_tracker,pa.t_in_fract)){return false;}

		// 2. Updating of particle characteristics
		pa.M = M_out;
		if (particle_time==-1){pa.t=-1;}
		else{
			pa.t += particle_time;
			pa.t_in_fract+=particle_time;
			pa.L_in_fract+=fabs(current_mesh.velocity)*advection_time;
		}
		pair<pointcpp<double>,pointcpp<double> > seg_current=make_pair(M_init, M_out);
		Segment_Particle_Time[seg_current].insert(pa.t); // when the particle reaches another position
	 }
	return true;
}
// Modified by Wenyu on 2026/1/8: Infinite matrix displacement with distance tracking
bool Transport::infinite_matrix_displacement_step(Particle & pa, double t_target, std::map<int,double>& moved_distances){
	// 0. Variables
	FractureMesh current_mesh = net_mesh.return_mesh(pa.mesh_index);
	static std::set<int> warned_zero_velocity;
	if ((current_mesh.velocity==0.0 || current_mesh.aperture==0.0) &&
	    warned_zero_velocity.insert(pa.mesh_index).second){
		cout << "Attention: in infinite_matrix_displacement_step aperture closed == 0: mesh_index="
		     << pa.mesh_index << " velocity=" << current_mesh.velocity
		     << " aperture=" << current_mesh.aperture << endl;
	}
	// If velocity is zero, just advance time to t_target; not stop simulation
	if (fabs(current_mesh.velocity== 0.0)){
		pointcpp<double> M_end(current_mesh.p_tar.p);
		pointcpp<double> M_start(current_mesh.p_ori.p);
		if (!HasConnectedMeshAtNode(net_mesh,pa.mesh_index,current_mesh.p_tar.index)){
			pa.M = M_start;
		}
		else{
			pa.M = M_end;
		}  
		double dt_remain = t_target - pa.t;
		if (dt_remain>0.0){pa.t_in_fract+=dt_remain;}
		pa.t = t_target;
		return true;
	}
	
	if (current_mesh.velocity<0&&pa.M.identic(pointcpp<double>(current_mesh.p_ori.p),EPSILON)){
		cout << "WARNING in infinite_matrix_displacement_step: flow velocity negative from origin to target" << endl;return false;
	}
	
	pointcpp<double> M_init = pa.M, M_out;
	current_mesh.distance_from_M_to_extremity(M_init, M_out);
	double advection_time = current_mesh.Mesh_Time_Advection(M_init, M_out);
	double particle_time;

	if (pa.t==-1){pa.M = M_out; return true;}
	if (!current_mesh.Get_Total_Time_From_Advec_Time1D(particle_time,advection_time,phys_param.Dm,phys_param.porosity,rng_tracker,pa.t_in_fract)){return false;}

	if (pa.t + particle_time > t_target){
		double dt_remain = t_target - pa.t;
		double adv_time_used = std::min(dt_remain, advection_time);
		pointcpp<double> M_new = current_mesh.Mesh_Scale_Advection(M_init, adv_time_used);
		pa.M = M_new;
		pa.t = t_target;
		pa.t_in_fract+=dt_remain;
		pa.L_in_fract+=fabs(current_mesh.velocity)*adv_time_used;
		moved_distances[pa.mesh_index]+=fabs(current_mesh.velocity)*adv_time_used;
		return true;
	}

	pa.M = M_out;
	if (particle_time==-1){pa.t=-1;}
	else{
		pa.t += particle_time;
		pa.t_in_fract+=particle_time;
		pa.L_in_fract+=fabs(current_mesh.velocity)*advection_time;
		moved_distances[pa.mesh_index]+=fabs(current_mesh.velocity)*advection_time;
	}
	pair<pointcpp<double>,pointcpp<double> > seg_current=make_pair(M_init, M_out);
	Segment_Particle_Time[seg_current].insert(pa.t); // when the particle reaches another position
	return true;
}

// ADDITION OF NEIGHBOURS EVEN IF THEIR FLUX IS LOWER THAN EPS_POINT_OUT_LIMIT_NO_FLOW TO OBSERVE DIFFUSION IN FRACTURED MEDIA (DR)
// affect the probabilities depending on the possible path with respect to the streamline routine method
bool Transport::select_neighbour_slr(vector<int>& possible_neighbours,vector<double>& possible_probabilities,FractureMesh Mesh_current,CgalPoint2D orip){
	size_t nb_pos = 0;				// Number of possible pathways from "ori" (only segments in which flux leaves "ori")
	double total_flow = 0.;			// Total flow leaving "ori"
	// the possible probabilities are filled with flow values, so compute the total flow and count the number of effectively possible neighbours (flow > 0)
	for(size_t i=0 ; i<possible_probabilities.size() ; i++){
		if(possible_probabilities[i]>0){		//#JRTRACKER+ROMAIN:PRECISION en train de retirer les bras morts (correction pr�alable des flux pour avoir un flux nul sur les ar�tes "bras morts")
			nb_pos++;
			total_flow += possible_probabilities[i];
		}
	}
	// SLR as a function of the number of possible paths
	switch(nb_pos){
		case(0) : { // limit point, the particle must go outside the system
			if(!this->domain.on_output_limit(orip)){
				cout << "WARNINIG in select_neighbour_slr (Transport.cpp): a point with no fluxes leaving the segment edge ori" << endl;
				return false;
			}
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				possible_probabilities[i] = 1;
			}
			return true;
				  } // case 0
		case(1) : { // only one possible path, its probability is set to one
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				if(possible_probabilities[i]>0)
					possible_probabilities[i] = 1;
				else
					possible_probabilities[i] = 0;
			}
				return true;
				  } // case 1
		case(2) : { // two possibilities : continue in the same fracture or turn (case 1) or turn left or right (case 2)
			// check if we are in case 1 or 2
			int q1 = -1; int q4 = -1;
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				CgalVector2D vect1(Mesh_current.p_ori.p, Mesh_current.p_tar.p);
				int index_neigh = possible_neighbours[i];
				FractureMesh mesh_neigh = net_mesh[index_neigh];
				CgalVector2D vect2(orip,mesh_neigh.p_tar.p);
				if((fabs((double)(vect1.x()/vect2.x()-vect1.y()/vect2.y()))<1e-10) && (possible_probabilities[i]>0))
					q1 = i; // index of the continuation in the same fracture
				else if(possible_probabilities[i]>0)
					q4 = i; // index of the bifurcation
			}
			if(q1!=-1 && q4 != -1){ // we find a continuation and a bifurcation, so apply SLR
				if(possible_probabilities[q1]<possible_probabilities[q4]){	// bifurcation dominates, so goes into it
					for(size_t i=0 ; i<possible_probabilities.size() ; i++)
						possible_probabilities[i] = 0.;
					possible_probabilities[q4] = 1;
				}else{														// straight line dominate, mixing probabilities
					double pq4 = possible_probabilities[q4]/possible_probabilities[q1];
					double pq1 = 1.-pq4;
					for(size_t i=0 ; i<possible_probabilities.size() ; i++)
						possible_probabilities[i] = 0.;
					possible_probabilities[q1] = pq1;
					possible_probabilities[q4] = pq4;
				}
			}else{ // we did not find a straight line, so there is two bifurcations, flow-equal probabilities
				for(size_t i=0 ; i<possible_probabilities.size() ; i++){
					if(possible_probabilities[i]>0)
						possible_probabilities[i] /= total_flow;
					else
						possible_probabilities[i] = 0;
				}
			}
				return true;
				  } // case 2
		case(3) : { // all path are possible excepting going back: probabilities = prorata of fluxes
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				if(possible_probabilities[i]>0)
					possible_probabilities[i] /= total_flow;
				else
					possible_probabilities[i] = 0;
			}
				return true;
				  } // case 3
		default : { // more than 3 possibilities, use prorata of fluxes
			for(size_t i=0 ; i<possible_probabilities.size() ; i++){
				if(possible_probabilities[i]>0)
					possible_probabilities[i] /= total_flow;
				else
					possible_probabilities[i] = 0;
			}
			return true;
				  }
	}
	return false;
}

// return the index of the next mesh from the current position and the current mesh of the particle
// modified to return false if the neighbor has been found
//int Transport::Mesh_Neighbour_Choice_And_Test(pointcpp<double> particle_pos, FractureMesh Mesh_current, RngStream_a& rng){
bool Transport::Mesh_Neighbour_Choice_And_Test(pointcpp<double> particle_pos,FractureMesh Mesh_current,
                                               RngStream_a& rng,int& neigh_index,bool& no_downstream){

	// 1. Variables
	CgalPoint2D current_pos(particle_pos.i,particle_pos.j);
	no_downstream = false;
	// 2. Determination of the corresponding mesh extremity
	FluxPoint2D ori;
	if(CGAL::squared_distance(Mesh_current.p_ori.p,current_pos)<CGAL::squared_distance(Mesh_current.p_tar.p,current_pos)){ori = Mesh_current.p_ori;}
	else{ori = Mesh_current.p_tar;}	// Transforms Flux_point_2D_neighbour in Flux_point_2D
	// 3. Determination of neighboring meshes and the associated probabilities
	vector<int> possible_neighbours;
	vector<double> possible_probabilities;
	vector<int> neighbours = Mesh_current.neigh_meshes;
	// Sorts among the neighbors of "Mesh_current->p_tar" and "Mesh_current->p_ori" the neighbors of "ori"
	size_t neighbours_size=neighbours.size() ;
	for(size_t i=0 ; i<neighbours_size ; i++){
		FractureMesh mesh_neigh = net_mesh[neighbours[i]];
		FluxPoint2D origin = mesh_neigh.p_ori;
		FluxPoint2D target = mesh_neigh.p_tar;
		if (origin.index == ori.index){
			possible_neighbours.push_back(neighbours[i]);
			double flux = net_mesh[neighbours[i]].velocity*net_mesh[neighbours[i]].aperture;		// Gets the flow of the mesh that will be used do determine probabilities
			possible_probabilities.push_back(flux);
		}
		else if (target.index == ori.index){
			possible_neighbours.push_back(neighbours[i]);
			double flux = -net_mesh[neighbours[i]].velocity*net_mesh[neighbours[i]].aperture;		// Gets the flow of the mesh that will be used do determine probabilities
			possible_probabilities.push_back(flux);
		}
	}
	// If no neighbors for the point
	if (possible_neighbours.empty()){
		std::cout << "WARNINIG in Mesh_Neighbour_Choice_And_Test (Transport.cpp): a point is connected to nothing" << endl;
		no_downstream = true;
		return false;
	}

	// 4. Computation of probabilities for the different possible outlet from the inlet point "ori"
	if (!select_neighbour_slr(possible_neighbours,possible_probabilities,Mesh_current,ori.p)){
		no_downstream = true;
		return false;
	}

	// 5. Choice of the outlet mesh from the probability distribution
	// Faire une fonction dans RNGstream_a.h
	// cumulative for possible probabilities
	int index = rng.drawPdf(possible_probabilities);
	neigh_index=possible_neighbours[index];
	return true;
	//return possible_neighbours[index];
	
}


// displacement until the right border of the domain
bool Transport::particle_displacement(Particle & pa){

	// 0. Variables
	int simu_option = this->num_param.simu_option;
	int new_mesh_index = -1; bool cont_disp=true;
	int nb_iter=0,nb_max_iter=1e4;

	// cout << "Initial position of the particle: "; print(pa.M.CgalPoint());

	// 1. Particle displacement until the right boundary of the domain
	while (!domain.on_output_limit(pa.M)&&cont_disp&&domain.IsInDomain(pa.M.CgalPoint())&&nb_iter<nb_max_iter){	// MODIF DR - 2020/11/03
		nb_iter++;
		// 1.1. Local displacement until a segment extremity
		// 1.1.1. Infinite matrix case
		if (simu_option==INFINITE_MATRIX){cont_disp=infinite_matrix_displacement(pa);}
		// 1.1.2. Finite matrix case
		else if (simu_option==FINITE_MATRIX){cont_disp=finite_matrix_displacement(pa);}
		else{cout << "WARNING in particle_displacement (Transport.cpp): simulation option not defined" << endl;}
		// 1.2. Choice of the next fracture
		// particle did not reach the end of the domain -> determination of the next mesh
		if (!domain.on_output_limit(pa.M)&&cont_disp&&domain.IsInDomain(pa.M.CgalPoint())){	// MODIF DR - 2020/11/03
			// determination of the next fracture from proportional flow
			FractureMesh mesh_current = net_mesh[pa.mesh_index];
			
			// modified by DR on 2025/12/11: before entering into a new fracture mesh, the aperture of the current fracture mesh is modified in relation with the time spent by the particle in this mesh
			 /*net_mesh_modified.ChangeAperture(pa.mesh_index,pa.t_in_fract,num_param.nb_part);
			 int Nt = num_param.nb_part;*/
			
			// modified by Wenyu on 2025/12/30: Count particles passed into fracture mesh and update the mesh grids
			// update aperture using empirical law
			net_mesh_modified.Nb_passed[pa.mesh_index]++;
			double Nb = (double)net_mesh_modified.Nb_passed[pa.mesh_index];
			double Nt = (double)num_param.nb_part;

			// --- update aperture ---
			net_mesh_modified.ChangeAperture(pa.mesh_index, pa.t_in_fract, Nb, Nt);
			
			
			bool no_downstream = false;
			if (!Mesh_Neighbour_Choice_And_Test(pa.M,mesh_current,rng_tracker,new_mesh_index,no_downstream)){
				if (no_downstream){
					std::set<int> connected_nodes = net_mesh.return_connected_nodes();
					ReassignToNearestConnectedEndpoint(pa,net_mesh,domain,connected_nodes,pa.t,false);
					return true;
				}
				return false;
			}

			// cout << "Particle moves from mesh " << pa.mesh_index << " to mesh " << new_mesh_index << endl;

			// particle updating
			pa.mesh_index = new_mesh_index;
			// when a particle enters into a new fracture, L_in_fract is set to the distance between the current position of the particle and the origin of the fracture
			if (mesh_current.fracture_index!=net_mesh[new_mesh_index].fracture_index){
				pa.L_in_fract=pa.M.Point_Distance(net_mesh.fractures[net_mesh[new_mesh_index].fracture_index].first.p);
			}
		}
	}

	// modified by DR on 2025/12/11: before entering into a new fracture mesh, the aperture of the current fracture mesh is modified in relation with the time spent by the particle in this mesh
	// modified by Wenyu on 2025/12/30.
	net_mesh_modified.Nb_passed[pa.mesh_index]++;
	double Nb = (double)net_mesh_modified.Nb_passed[pa.mesh_index];
	double Nt = (double)num_param.nb_part;

	net_mesh_modified.ChangeAperture(pa.mesh_index, pa.t_in_fract, Nb, Nt);
        if (nb_iter==nb_max_iter){cout << "Particle stopped because of nb_max_iter, simulation stopped"; pa.t=-1; return false;}

	return cont_disp;
}
// Modified by Wenyu on 2026/1/8: Particle displacement with distance tracking
bool Transport::particle_displacement_step(Particle & pa, double dt_step, std::map<int,double>& moved_distances){
	double t_target = pa.t + dt_step;
	int simu_option = this->num_param.simu_option;
	int new_mesh_index = -1; bool cont_disp=true;
	int nb_iter=0,nb_max_iter=1e4;

	while (!domain.on_output_limit(pa.M)&&cont_disp&&domain.IsInDomain(pa.M.CgalPoint())&&pa.t<t_target&&nb_iter<nb_max_iter){
		nb_iter++;
		if (simu_option==INFINITE_MATRIX){cont_disp=infinite_matrix_displacement_step(pa,t_target,moved_distances);}
		else if (simu_option==FINITE_MATRIX){cont_disp=finite_matrix_displacement_step(pa,t_target,moved_distances);}
		else{cout << "WARNING in particle_displacement_step (Transport.cpp): simulation option not defined" << endl;}
		if (pa.t>=t_target){break;}
		if (!domain.on_output_limit(pa.M)&&cont_disp&&domain.IsInDomain(pa.M.CgalPoint())){
			FractureMesh mesh_current = net_mesh[pa.mesh_index];
			bool no_downstream = false;
			if (!Mesh_Neighbour_Choice_And_Test(pa.M,mesh_current,rng_tracker,new_mesh_index,no_downstream)){
				if (no_downstream){
					std::set<int> connected_nodes = net_mesh.return_connected_nodes();
					ReassignToNearestConnectedEndpoint(pa,net_mesh,domain,connected_nodes,t_target,true);
					break;
				}
				return false;
			}
			pa.mesh_history.push_back(pa.mesh_index);
			pa.prev_mesh_index = pa.mesh_index;
			pa.mesh_index = new_mesh_index;
			if (mesh_current.fracture_index!=net_mesh[new_mesh_index].fracture_index){
				pa.L_in_fract=pa.M.Point_Distance(net_mesh.fractures[net_mesh[new_mesh_index].fracture_index].first.p);
			}
		}
	}

	if (nb_iter==nb_max_iter){cout << "Particle stopped because of nb_max_iter, simulation stopped"; pa.t=-1; return false;}
	return cont_disp;
}

bool Transport::Particles_Transport(map<int,double> & arrival_times,int option_injection){
	// 0. Variables
	int nb_part=num_param.nb_part;
	// 1. Particles initialization and injection on the left border
	vector<Particle> part_vect = Particles_Injection(option_injection);
	// 2. Particle displacement
	// Modified by Wenyu on 2026/1/8: use reaction time step and output interval
	double dt_step = num_param.reaction_dt;
	if (dt_step<=0.0){
		cout << "WARNING in Particles_Transport (Transport.cpp): invalid reaction time step" << endl;
		return false;
	}
	double output_interval = num_param.output_interval;
	if (output_interval<=0.0){output_interval = dt_step;}
	double reference_injected_count_dt = 0.0;
	if (num_param.t_injection>0.0){
		if (nb_part>1){
			reference_injected_count_dt = dt_step * ((double)(nb_part-1) / num_param.t_injection);
		}
		else if (nb_part==1){
			reference_injected_count_dt = dt_step / num_param.t_injection;
		}
	}
	else{
		// Instantaneous injection: use the total injected particles as the reference scale.
		reference_injected_count_dt = (double)nb_part;
	}
	if (reference_injected_count_dt<=0.0){
		cout << "WARNING in Particles_Transport (Transport.cpp): invalid reference injected count per reaction step" << endl;
		return false;
	}

	double t_current = 0.0;
	double next_output_time = 0.0;
	RecordPositions(0.0,part_vect);
	WriteDFNSnapshot(net_mesh_modified,param_full.code_path,0,0.0);
	next_output_time += output_interval;
	int max_output_steps = (output_interval>0.0)? (int)floor(num_param.t_injection/output_interval)+1 : 0;
	int recorded_steps = 1;

	int step_index = 0;
	int total_steps_est = (dt_step>0.0)? (int)floor(num_param.t_max/dt_step)+1 : 0;
	int timing_prints = 0;
	int timing_total = (total_steps_est>0)? std::min(10,total_steps_est) : 0;
	double total_step_time = 0.0;
	std::map<int,double> particle_pass_count;
	std::set<int> zero_aperture_logged;
	// Modified by Wenyu on 2026/1/13: loop over reaction time steps
	{
		std::map<int,pointcpp<double> > init_pos;
		CollectInputPositions(net_mesh_modified,domain,init_pos);
		if (option_injection==1 && !init_pos.empty()){
			double min_dist = 1e9;
			pointcpp<double> middle_left_border = domain.return_middle_input_limit();
			int init_mesh_index = -1;
			pointcpp<double> init_pt_single;
			for (map<int,pointcpp<double> >::iterator it=init_pos.begin();it!=init_pos.end();it++){
				double dist = middle_left_border.Point_Distance(it->second);
				if (dist<min_dist){
					min_dist = dist;
					init_mesh_index = it->first;
					init_pt_single = it->second;
				}
			}
			init_pos.clear();
			if (init_mesh_index!=-1){
				init_pos[init_mesh_index] = init_pt_single;
			}
		}
		std::vector<int> mesh_ids;
		std::vector<pointcpp<double> > positions;
		for (map<int,pointcpp<double> >::iterator it=init_pos.begin();it!=init_pos.end();it++){
			mesh_ids.push_back(it->first);
			positions.push_back(it->second);
		}
		std::vector<int> counts;
		if (!mesh_ids.empty()){
			if (option_injection==2){
				std::vector<double> weights(mesh_ids.size(),0.0);
				for (size_t i=0;i<mesh_ids.size();i++){
					FractureMesh mesh = net_mesh_modified.return_mesh(mesh_ids[i]);
					double w = mesh.aperture*std::fabs(mesh.velocity);
					weights[i] = (w>0.0)? w : 0.0;
				}
				if (ComputeInjectionCounts(nb_part,weights,counts)){
					PrintInjectionCounts("Initial inlet allocation",mesh_ids,counts);
				}
			}
			else{
				std::vector<double> weights(mesh_ids.size(),1.0);
				if (ComputeInjectionCounts(nb_part,weights,counts)){
					PrintInjectionCounts("Initial inlet allocation",mesh_ids,counts);
				}
			}
		}
	}
	while (t_current < num_param.t_max+EPSILON){
		auto step_start = std::chrono::steady_clock::now();
		bool any_active = false;
		bool rebuild_needed = false;
		std::map<int,int> step_full_count;
		std::map<int,double> step_partial_sum;
		std::vector<int> ready_indices;
		for (int i=0;i<nb_part;i++){
			if (part_vect[i].t==-1){continue;}
			if (part_vect[i].mesh_index==-1 && part_vect[i].t<=t_current+EPSILON){
				ready_indices.push_back(i);
			}
		}
		if (!ready_indices.empty()){
			std::map<int,pointcpp<double> > init_pos;
			CollectInputPositions(net_mesh_modified,domain,init_pos);
			if (option_injection==1 && !init_pos.empty()){
				double min_dist = 1e9;
				pointcpp<double> middle_left_border = domain.return_middle_input_limit();
				int init_mesh_index = -1;
				pointcpp<double> init_pt_single;
				for (map<int,pointcpp<double> >::iterator it=init_pos.begin();it!=init_pos.end();it++){
					double dist = middle_left_border.Point_Distance(it->second);
					if (dist<min_dist){
						min_dist = dist;
						init_mesh_index = it->first;
						init_pt_single = it->second;
					}
				}
				init_pos.clear();
				if (init_mesh_index!=-1){
					init_pos[init_mesh_index] = init_pt_single;
				}
			}
			int nb_pos = init_pos.size();
			if (nb_pos==0){
				cout << "WARNING in Particles_Transport (Transport.cpp): no active inlet mesh for injection at step="
				     << (step_index+1) << " t=" << t_current << endl;
				return false;
			}
			else{
				std::vector<int> mesh_ids;
				std::vector<pointcpp<double> > positions;
				mesh_ids.reserve(nb_pos);
				positions.reserve(nb_pos);
				for (map<int,pointcpp<double> >::iterator it=init_pos.begin();it!=init_pos.end();it++){
					mesh_ids.push_back(it->first);
					positions.push_back(it->second);
				}
				std::vector<int> counts;
				if (option_injection==2){
					std::vector<double> weights(nb_pos,0.0);
					for (int i=0;i<nb_pos;i++){
						FractureMesh mesh = net_mesh_modified.return_mesh(mesh_ids[i]);
						double w = mesh.aperture*std::fabs(mesh.velocity);
						weights[i] = (w>0.0)? w : 0.0;
					}
					if (!ComputeInjectionCounts((int)ready_indices.size(),weights,counts)){
						cout << "WARNING in Particles_Transport (Transport.cpp): no positive inlet flux for injection at step="
						     << (step_index+1) << " t=" << t_current << endl;
						cout << "Inlet mesh aperture/velocity:" << endl;
						for (int i=0;i<nb_pos;i++){
							FractureMesh mesh = net_mesh_modified.return_mesh(mesh_ids[i]);
							cout << "mesh_index=" << mesh_ids[i]
							     << " aperture=" << mesh.aperture
							     << " velocity=" << mesh.velocity << endl;
						}
						return false;
					}
				}
				else{
					std::vector<double> weights(nb_pos,1.0);
					ComputeInjectionCounts((int)ready_indices.size(),weights,counts);
				}
				if (nb_pos>0){
					int idx = 0;
					for (int i=0;i<nb_pos;i++){
			for (int c=0;c<counts[i] && idx<(int)ready_indices.size();c++){
				Particle & pa = part_vect[ready_indices[idx]];
				pa.M = positions[i];
				pa.mesh_index = mesh_ids[i];
				pa.t_injection = pa.t;
							pa.prev_mesh_index = -1;
							pa.mesh_history.clear();
							pa.intersection_history.clear();
							pa.L_in_fract = 0.0;
							pa.t_in_fract = 0.0;
							pa.t_in_fract_prev = 0.0;
							Initial_Position_Particles[pa.no] = pa.M;
							idx++;
						}
					}
				}
			}
		}
		for (int i=0;i<nb_part;i++){
			Particle & pa = part_vect[i];
			if (!std::isfinite(pa.t)){
				cout << "WARNING in Particles_Transport (Transport.cpp): particle time is NaN, reassign to inlet: particle="
				     << pa.no << " mesh_index=" << pa.mesh_index << endl;
				ResetParticleForInlet(pa,t_current);
				continue;
			}
			if (pa.t==-1){continue;}
			if (pa.mesh_index==-1){continue;}
			if (pa.t>t_current+EPSILON){continue;}
			any_active = true;
			std::map<int,double> moved_distances;
			if (!particle_displacement_step(pa,dt_step,moved_distances)){return false;}
			if (pa.t!=-1 && domain.on_output_limit(pa.M)){
				if (std::isfinite(pa.t)){
					arrival_times[pa.no]=pa.t;
					pa.t=-1;
				}
				else{
					cout << "WARNING in Particles_Transport (Transport.cpp): particle time NaN at output, reassign to inlet: particle="
					     << pa.no << " mesh_index=" << pa.mesh_index << endl;
					ResetParticleForInlet(pa,t_current);
				}
			}
			for (std::map<int,double>::iterator it= moved_distances.begin(); it!=moved_distances.end(); it++){
				int mesh_id = it->first;
				double length = net_mesh_modified.return_mesh(mesh_id).ReturnLength();
				double dist = it->second;
				if (!std::isfinite(length) || length<=0.0){
					cout << "WARNING in Particles_Transport (Transport.cpp): invalid mesh length="
					     << length << " mesh_index=" << mesh_id << endl;
					continue;
				}
				if (!std::isfinite(dist)){
					cout << "WARNING in Particles_Transport (Transport.cpp): invalid moved distance="
					     << dist << " mesh_index=" << mesh_id << endl;
					continue;
				}
				double fraction = dist/length;
				if (!std::isfinite(fraction)){
					cout << "WARNING in Particles_Transport (Transport.cpp): invalid fraction="
					     << fraction << " mesh_index=" << mesh_id << endl;
					continue;
				}
				int full = (int)floor(fraction+EPSILON);
				if (full>=1){
					step_full_count[mesh_id] += full;
				}
				else{
					step_full_count[mesh_id] += 0;
					step_partial_sum[mesh_id] += fraction;
				}
			}
		}

		for (std::map<int,int>::iterator it= step_full_count.begin(); it!=step_full_count.end(); it++){
			int mesh_id = it->first;
			double Nb_eff = (double)step_full_count[mesh_id] + step_partial_sum[mesh_id];
			if (!std::isfinite(Nb_eff)){
				cout << "WARNING in Particles_Transport (Transport.cpp): Nb_eff is NaN/inf for mesh_index="
				     << mesh_id << endl;
				continue;
			}
			particle_pass_count[mesh_id] += Nb_eff;
			double b_old = net_mesh_modified.return_mesh(mesh_id).aperture;
			net_mesh_modified.ChangeAperture(mesh_id,dt_step,Nb_eff,reference_injected_count_dt);
			double b_new = net_mesh_modified.return_mesh(mesh_id).aperture;
			if (b_old>0.0 && b_new==0.0){
				rebuild_needed = true;
			}
			if (b_old>0.0 && b_new==0.0 && zero_aperture_logged.insert(mesh_id).second){
				cout << "WARNING: aperture->0 at step=" << (step_index+1)
				     << " t=" << (t_current+dt_step)
				     << ", mesh_index=" << mesh_id
				     << ", aperture_old=" << b_old
				     << ", particles_passed=" << particle_pass_count[mesh_id] << endl;
			}
		}

		t_current += dt_step;
		step_index++;
		if (t_current+EPSILON >= next_output_time && recorded_steps<max_output_steps){
			RecordPositions(t_current,part_vect);
			WriteDFNSnapshot(net_mesh_modified,param_full.code_path,recorded_steps,t_current);
			next_output_time += output_interval;
			recorded_steps++;
		}
		// If rebuild is needed, recalculate the conectivity and flow field
		if (!rebuild_needed){
			for (size_t i=0;i<net_mesh_modified.meshes.size();i++){
				if (net_mesh_modified.meshes[i].aperture<=0.0){
					rebuild_needed = true;
					break;
				}
			}
		}
		// Rebuild the network if some apertures have reached zero（from variable not recalculate by function）
		if (rebuild_needed){
			NetworkMeshes mesh_before = net_mesh_modified;
			net_mesh_before_backbone = mesh_before;
			std::map<int,int> old_to_new;
			std::vector<int> removed_meshes;
			double rebuild_seconds = 0.0;
			std::map<int,int> node_old_to_new;
			NetworkMeshes tmp_mesh = net_mesh_modified;
			RebuildNetworkRemovingZeroAperture(tmp_mesh,old_to_new,node_old_to_new,removed_meshes,rebuild_seconds);
			cout << "Timing: rebuild network after aperture=0 took " << rebuild_seconds << " s" << endl;
			// Rebuild backbone on the updated network
			NetworkMeshes backbone_mesh = tmp_mesh.return_backbone(param_full,EPSILON);
			if (backbone_mesh.meshes.size()==0){
				cout << "WARNING: DFN disconnected after rebuild at step=" << (step_index+1)
				     << " t=" << (t_current+dt_step) << endl;
				return false;
			}
			FractureMeshNumbering(backbone_mesh);
			UpdateCptFract(backbone_mesh);
			net_mesh_modified = backbone_mesh;
			std::map<int,int> backbone_old_to_new;
			RebuildNodeNumbering(net_mesh_modified,backbone_old_to_new);
			std::set<int> connected_nodes = net_mesh_modified.return_connected_nodes();
			for (size_t i=0;i<part_vect.size();i++){
				if (part_vect[i].t==-1){continue;}
				if (part_vect[i].mesh_index==-1){continue;}
				FractureMesh old_mesh = mesh_before.return_mesh(part_vect[i].mesh_index);
				int new_index = FindMatchingMeshIndex(net_mesh_modified,old_mesh);
				if (new_index<0){
					RemoveParticleFromSystem(part_vect[i],(t_current+dt_step),"mesh_removed");
					continue;
				}
				part_vect[i].mesh_index = new_index;
				part_vect[i].mesh_history.clear();
				part_vect[i].prev_mesh_index = -1;
				FractureMesh cur_mesh = net_mesh_modified.return_mesh(part_vect[i].mesh_index);
				if (connected_nodes.find(cur_mesh.p_ori.index)==connected_nodes.end() ||
				    connected_nodes.find(cur_mesh.p_tar.index)==connected_nodes.end()){
					RemoveParticleFromSystem(part_vect[i],(t_current+dt_step),"mesh_disconnected");
				}
			}
		}
		if (any_active){
			ComputeFlowVelocities(net_mesh_modified,param_full,domain);
			NormalizeMeshDirections(net_mesh_modified);
			if (ZeroApertureIfLowVelocity(net_mesh_modified)){
				rebuild_needed = true;
			}
			net_mesh = net_mesh_modified;
			auto step_end = std::chrono::steady_clock::now();
			double step_sec = std::chrono::duration<double>(step_end-step_start).count();
			total_step_time += step_sec;
			double min_aperture = 0.0;
			int min_aperture_mesh = -1;
			if (!net_mesh_modified.meshes.empty()){
				min_aperture = net_mesh_modified.meshes[0].aperture;
				min_aperture_mesh = net_mesh_modified.meshes[0].mesh_index;
				for (size_t i=1;i<net_mesh_modified.meshes.size();i++){
					if (net_mesh_modified.meshes[i].aperture<min_aperture){
						min_aperture = net_mesh_modified.meshes[i].aperture;
						min_aperture_mesh = net_mesh_modified.meshes[i].mesh_index;
					}
				}
			}
			if (timing_total>0){
				int target_index = (int)floor(((double)(timing_prints+1) * (double)total_steps_est) / (double)timing_total);
				if (step_index>=target_index){
					cout << "Timing step " << step_index << "/" << total_steps_est
					     << " t=" << t_current
					     << ": total_step=" << step_sec
					     << " s, min_aperture=" << min_aperture
					     << " mesh_index=" << min_aperture_mesh << endl;
					timing_prints++;
				}
			}
		}
		else{
			break;
		}
	}
	{
		std::map<int,pointcpp<double> > init_pos;
		CollectInputPositions(net_mesh_modified,domain,init_pos);
		if (option_injection==1 && !init_pos.empty()){
			double min_dist = 1e9;
			pointcpp<double> middle_left_border = domain.return_middle_input_limit();
			int init_mesh_index = -1;
			pointcpp<double> init_pt_single;
			for (map<int,pointcpp<double> >::iterator it=init_pos.begin();it!=init_pos.end();it++){
				double dist = middle_left_border.Point_Distance(it->second);
				if (dist<min_dist){
					min_dist = dist;
					init_mesh_index = it->first;
					init_pt_single = it->second;
				}
			}
			init_pos.clear();
			if (init_mesh_index!=-1){
				init_pos[init_mesh_index] = init_pt_single;
			}
		}
		std::vector<int> mesh_ids;
		std::vector<pointcpp<double> > positions;
		for (map<int,pointcpp<double> >::iterator it=init_pos.begin();it!=init_pos.end();it++){
			mesh_ids.push_back(it->first);
			positions.push_back(it->second);
		}
		std::vector<int> counts;
		if (!mesh_ids.empty()){
			if (option_injection==2){
				std::vector<double> weights(mesh_ids.size(),0.0);
				for (size_t i=0;i<mesh_ids.size();i++){
					FractureMesh mesh = net_mesh_modified.return_mesh(mesh_ids[i]);
					double w = mesh.aperture*std::fabs(mesh.velocity);
					weights[i] = (w>0.0)? w : 0.0;
				}
				if (ComputeInjectionCounts(nb_part,weights,counts)){
					PrintInjectionCounts("Final inlet allocation",mesh_ids,counts);
				}
			}
			else{
				std::vector<double> weights(mesh_ids.size(),1.0);
				if (ComputeInjectionCounts(nb_part,weights,counts)){
					PrintInjectionCounts("Final inlet allocation",mesh_ids,counts);
				}
			}
		}
	}
	cout << "Timing total for all steps: " << total_step_time << " s" << endl;
	return true;
}
// Modified by Wenyu on 2026/1/8: Record particle positions at specific time
void Transport::RecordPositions(double time, const std::vector<Particle>& parts){
	std::vector<ParticleSnapshot> snapshots;
	snapshots.reserve(parts.size());
	for (size_t i=0;i<parts.size();i++){
		if (parts[i].t==-1){continue;}
		if (parts[i].mesh_index==-1){continue;}
		if (parts[i].t>time+EPSILON){continue;}
		ParticleSnapshot snap;
		snap.no = parts[i].no;
		snap.mesh_index = parts[i].mesh_index;
	if (parts[i].mesh_index>=0){
		snap.mesh_index_original = net_mesh_modified.return_mesh(parts[i].mesh_index).original_mesh_index;
	}
		else{
			snap.mesh_index_original = -1;
		}
		snap.t_injection = parts[i].t_injection;
		snap.M = parts[i].M;
		snapshots.push_back(snap);
	}
	Position_Time_Particles[time]=snapshots;
}
// Modified by Wenyu on 2026/1/8: Write particle position snapshots to CSV files
void Transport::WritePositionSnapshotsCSV(const std::string& output_path) const{
	int file_index = 0;
	for (std::map<double,std::vector<ParticleSnapshot> >::const_iterator it=Position_Time_Particles.begin();
	it!=Position_Time_Particles.end(); it++){
		std::ostringstream name_stream;
		name_stream << output_path << "/Output/particle_positions_t"
		            << std::setw(6) << std::setfill('0') << file_index << ".csv";
		std::string file_name = name_stream.str();
		std::ofstream output(file_name.c_str(), std::ofstream::out);
		if (!output.is_open()){
			cout << "WARNING in WritePositionSnapshotsCSV (Transport.cpp): cannot open file " << file_name << endl;
			return;
		}
		double t = it->first;
		const std::vector<ParticleSnapshot>& snaps = it->second;
		output << "time,particle_id,x,y,mesh_index,mesh_index_original" << endl;
		for (size_t i=0;i<snaps.size();i++){
			output << t << "," << snaps[i].no << "," << snaps[i].M.i << "," << snaps[i].M.j
			       << "," << snaps[i].mesh_index
			       << "," << snaps[i].mesh_index_original << endl;
		}
		output.close();
		file_index++;
	}
}

void Transport::PrintTravelTimeSamplesAtOutputIntervals(int max_samples) const{
	for (std::map<double,std::vector<ParticleSnapshot> >::const_iterator it=Position_Time_Particles.begin();
	     it!=Position_Time_Particles.end(); it++){
		double t = it->first;
		const std::vector<ParticleSnapshot>& snaps = it->second;
		cout << "Diagnostics sample travel time at t=" << t << ":" << endl;
		int printed = 0;
		for (size_t i=0;i<snaps.size() && printed<max_samples;i++){
			double travel = t - snaps[i].t_injection;
			cout << "particle=" << snaps[i].no
			     << " travel_time=" << travel
			     << " mesh_index=" << snaps[i].mesh_index << endl;
			printed++;
		}
	}
}
