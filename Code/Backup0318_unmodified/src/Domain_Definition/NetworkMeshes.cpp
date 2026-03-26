/*
 * NetworkMeshes.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#include "NetworkMeshes.h"
#include "Domain.h"
#include "DFNComputation.h"
#include "../Utilitaries/Critical_Operations.h"
#include "../Utilitaries/Structures.h"
#include "../Utilitaries/UblasStructures.h"
#include "../Utilitaries/LinearSystem.h"
#include "../Chemistry/Chemistry.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <boost/math/special_functions/erf.hpp>

using namespace std;

void NodeInsertionMap(int,int,FractureMesh,NodeFracturesMap&);
void add_neigh_nodes(set<int>&,NodeFracturesMap,int);
// Modified by Wenyu on 2026/1/7; check if left and right borders are connected without going through the filled mesh
static bool LeftRightConnectedWithoutMesh(const NetworkMeshes& net_mesh,int skip_mesh_index){
	NodeFracturesMap node_fract_map;
	for (vector<FractureMesh>::const_iterator it=net_mesh.meshes.begin();it!=net_mesh.meshes.end();it++){
		if (it->mesh_index==skip_mesh_index){continue;}
		if (it->aperture<=0.0){continue;}
		NodeInsertionMap(it->p_ori.index,it->p_tar.index,*it,node_fract_map);
		NodeInsertionMap(it->p_tar.index,it->p_ori.index,*it,node_fract_map);
	}
	set<int> connected_nodes;
	for (BordersMap::const_iterator it=net_mesh.border_map.begin();it!=net_mesh.border_map.end();it++){
		if (it->second==LEFT_BORDER){
			connected_nodes.insert(it->first);
			add_neigh_nodes(connected_nodes,node_fract_map,it->first);
		}
	}
	if (connected_nodes.empty()){return false;}
	for (BordersMap::const_iterator it=net_mesh.border_map.begin();it!=net_mesh.border_map.end();it++){
		if (it->second==RIGHT_BORDER && connected_nodes.find(it->first)!=connected_nodes.end()){
			return true;
		}
	}
	return false;
}


NetworkMeshes::NetworkMeshes(){cpt_inter=-1;cpt_mesh=-1;cpt_fract=-1;}
NetworkMeshes::~NetworkMeshes(){}
NetworkMeshes::NetworkMeshes(string code_path,string file_name_DFN,Domain domain_){
	// 1. Variables
	domain=domain_;
	cpt_inter=-1;cpt_mesh=-1;cpt_fract=-1;
	string DFN_generation_option;
	// 2. Read the DFN generation option
	// Domain definition from files
	string file_name = code_path + "/Input/DFN_files/"+file_name_DFN;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){fichier >> DFN_generation_option;}
	fichier.close();
	// 3. Generate the fracture network
	if (DFN_generation_option=="file"){
		NetworkMeshes net_mesh=*this;
		NetworkMeshesFromFile(code_path,file_name_DFN,net_mesh);
		*this=net_mesh;
	}
	else if (DFN_generation_option=="generation_sierpinski"){
		NetworkMeshes net_mesh=*this;
		NetworkMeshesGenerationSierpinski(code_path,file_name_DFN,net_mesh);
		*this=net_mesh;
	}
	else if (DFN_generation_option=="generation_realistic"){
		NetworkMeshes net_mesh=*this;
		NetworkMeshesGenerationRealistic(code_path,file_name_DFN,net_mesh);
		*this=net_mesh;
	}
	else if (DFN_generation_option=="generation_realistic2"){
		NetworkMeshes net_mesh=*this;
		NetworkMeshesGenerationRealistic2(code_path,file_name_DFN,net_mesh);
		*this=net_mesh;
	}
	else if (DFN_generation_option=="generation_realistic3"){
                NetworkMeshes net_mesh=*this;
                NetworkMeshesGenerationRealistic3(code_path,file_name_DFN,net_mesh);
                *this=net_mesh;
        }
	else{cout << "WARNING in NetworkMeshes (NetworkMeshes.cpp): option not defined" << endl; cout << DFN_generation_option << endl;}
}

// fracture network definition for inversion study (large number of forward simulations)
//NetworkMeshes::NetworkMeshes(Domain domain_,double C_density,double D_dim,double coeff_theta1,double coeff_theta2,double r_min,double fract_aperture,int seed){
/*NetworkMeshes::NetworkMeshes(Domain domain_,double C_density,double D_dim,double r_min,double fract_aperture,int seed){
	// 1. Variables
	domain=domain_;
	cpt_inter=-1;cpt_mesh=-1;cpt_fract=-1;

	// 2. Fracture network properties
	// nb_fract is equal to the integer of i_fract rounded up.

	FractureMesh frac_mesh; Segment2D seg;double fract_length,angle;CgalPoint2D fract_center;
	double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();

	if (Lx!=Ly){cout << "PROBLEM in NetworkMeshes (NetworkMeshes.cpp): Lx and Ly must be equal" << endl; exit(0);}

	int nb_fract=std::ceil(C_density*std::pow(r_min/Lx,-(double)D_dim));
        int cpt_fract=0;

	NetworkMeshes net_mesh=*this;
	srand(seed);
	while (cpt_fract<nb_fract){
		cpt_fract++;
		// center of the fracture
		fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
		// length of the fracture
		fract_length=std::pow(C_density/(double)cpt_fract,1./(double)D_dim)*Lx;
		//angle of the fracture
		//if ((double)cpt_fract/2.==(int)cpt_fract/2){angle=coeff_theta1*PI/180;}
		//else {angle=coeff_theta2*PI/180;}
		angle = Uniform()*PI;
		// definition of a segment
		seg=Segment2D(fract_center,fract_length,angle);
		domain.SegmentIntersectDomain(seg,seg);
		frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
		AddFracture(frac_mesh,net_mesh);
	}
	// 4. Intersections determination
	ComputeIntersections(net_mesh);
	*this=net_mesh;
}*/

NetworkMeshes::NetworkMeshes(Domain domain_,double C_density,double D_dim,double r_min,double fract_aperture,int seed){
        // 1. Variables
        domain=domain_;
        cpt_inter=-1;cpt_mesh=-1;cpt_fract=-1;

        // 2. Fracture network properties
        // nb_fract is equal to the integer of i_fract rounded up.

        FractureMesh frac_mesh; Segment2D seg;double fract_length,angle;CgalPoint2D fract_center;
        double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();

        if (Lx!=Ly){cout << "PROBLEM in NetworkMeshes (NetworkMeshes.cpp): Lx and Ly must be equal" << endl; exit(0);}

        int nb_fract=std::ceil(C_density*std::pow(r_min/Lx,-(double)D_dim));
        int cpt_fract=0;

        NetworkMeshes net_mesh=*this;
        srand(seed);
        while (cpt_fract<nb_fract){
                cpt_fract++;
                // center of the fracture
                fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
                // length of the fracture
                fract_length=std::pow(C_density/(double)cpt_fract,1./(double)D_dim)*Lx;
                //angle of the fracture
                angle = Uniform()*PI;
                // definition of a segment
                seg=Segment2D(fract_center,fract_length,angle);
                domain.SegmentIntersectDomain(seg,seg);
                frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
                AddFracture(frac_mesh,net_mesh);
        }
        // 4. Intersections determination
        ComputeIntersections(net_mesh);
        *this=net_mesh;
}

// modified by DR on 2025/12/11: change the mesh aperture depending on the time spent by the particle in the mesh
// modified by Wenyu on 2025/12/30: Using new writtin function ComputeAperture to modified the aperture height of the mesh
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
void NetworkMeshes::ChangeAperture(int mesh_index, double dt_in_fract, double Nb, double Nt)
{
    if (mesh_index < 0 || mesh_index >= (int)meshes.size()) return;
    if (Nt <= 0) return;
    if (Nb < 0) Nb = 0;
    if (dt_in_fract < 0.0) dt_in_fract = 0.0;

    double b_old = meshes[mesh_index].aperture;
    double ratio = safe_ratio(Nb, Nt);
	double delta = K_REACTION * dt_in_fract * ratio;// k * t * (Nb/Nt), I didnt use the function in Chemistry.cpp to avoid extra function call overhead
    double b_new  = ComputeAperture(b_old, dt_in_fract, Nb, Nt, true); 

    if (b_new < 0.0) b_new = 0.0;

    // DEBUG BLOCK (Delete later)
    // ---------------------------
/*
    std::cerr << std::scientific << std::setprecision(6)
              << "[DBG] mesh=" << mesh_index
              << " b_old=" << b_old
              << " delta=" << delta
              << " b_new=" << b_new
              << " k=" << K_REACTION
              << " dt=" << dt_in_fract
              << " ratio=" << ratio
              << std::endl;*/
    // Write back
    // ---------------------------
    meshes[mesh_index].aperture = b_new;
	if (b_new==0.0){
		std::cout << "dt_in_fract=" << dt_in_fract << ", Nb=" << Nb << ", Nt=" << Nt << ", old aperture=" << b_old << ", new aperture=" << b_new << std::endl;
		std::cout << "ratio=" << ratio << std::endl;
		meshes[mesh_index].velocity = 0.0;
		if (!LeftRightConnectedWithoutMesh(*this,mesh_index)){
			cout << "WARNING in ChangeAperture (NetworkMeshes.cpp): aperture=0 may disconnect DFN; mesh_index="
			     << mesh_index << endl;
		}
	}
}
// modified by Wenyu on 2026/1/8: change the mesh aperture by a given ratio
void NetworkMeshes::ChangeApertureByRatio(int mesh_index, double travel_ratio)
{
	if (mesh_index < 0 || mesh_index >= (int)meshes.size()) return;
	if (travel_ratio <= 0.0) return;
	double b_old = meshes[mesh_index].aperture;
	double b_new = ComputeAperture(b_old, travel_ratio, 1, 1, true);
	meshes[mesh_index].aperture = b_new;
	if (b_new==0.0){
		meshes[mesh_index].velocity = 0.0;
		if (!LeftRightConnectedWithoutMesh(*this,mesh_index)){
			cout << "WARNING in ChangeApertureByRatio (NetworkMeshes.cpp): aperture=0 may disconnect DFN; mesh_index="
			     << mesh_index << endl;
		}
	}
}



int NetworkMeshes::ReturnIndex(CgalPoint2D pt){
	for (size_t i=0;i<meshes.size();i++){
		if (Points_Equal(pt,meshes[i].p_ori.p)){return meshes[i].p_ori.index;}
		if (Points_Equal(pt,meshes[i].p_tar.p)){return meshes[i].p_tar.index;}
	}
	return -1;
}

void Border_Definition(FluxPoint2D pt_inter,NetworkMeshes& net_mesh){
	string border=net_mesh.domain.ReturnBorder(pt_inter.p);
	if (border!=NO_BORDER){net_mesh.border_map[pt_inter.index]=border;}
	net_mesh.nodes_map[pt_inter.index]=pt_inter.p;
}

void AddFracture(FractureMesh frac_mesh,NetworkMeshes& net_mesh){
	// indices determination
	int cpt_index=net_mesh.cpt_inter;
	// definition of extremity numbering
	int cpt1=net_mesh.ReturnIndex(frac_mesh.p_ori.p),cpt2=net_mesh.ReturnIndex(frac_mesh.p_tar.p);
	if (cpt1==-1){cpt1=cpt_index+1;cpt_index++;}
	if (cpt2==-1){cpt2=cpt_index+1;cpt_index++;}
	FluxPoint2D pt1(frac_mesh.p_ori.p,cpt1),pt2(frac_mesh.p_tar.p,cpt2);
	//net_mesh.cpt_mesh=net_mesh.cpt_mesh+1;
	FractureMesh new_fract(pt1,pt2,frac_mesh.aperture,frac_mesh.fracture_index);//,net_mesh.cpt_mesh);
	net_mesh.meshes.push_back(new_fract);
	// boundary conditions update
	net_mesh.cpt_inter=cpt_index;
	Border_Definition(pt1,net_mesh);
	Border_Definition(pt2,net_mesh);
}

void DivideDomain(pointcpp<double> pt_min,pointcpp<double> pt_max,int nb_division,double aperture,NetworkMeshes& net_mesh,double Lx){
	// variables
	FractureMesh frac_mesh;int cpt_fract=net_mesh.cpt_fract;
	double x_current=pt_min.i,y_current=pt_min.j,x_min=pt_min.i,y_min=pt_min.j,x_max=pt_max.i,y_max=pt_max.j;
	double fract_spacing_x=(pt_max.i-pt_min.i)/nb_division,fract_spacing_y=(pt_max.j-pt_min.j)/nb_division;
	// creation of the horizontal fractures
	for (int i=0;i<nb_division-1;i++){
		y_current+=fract_spacing_y;
		cpt_fract++;
		frac_mesh=FractureMesh(FluxPoint2D(x_min,y_current),FluxPoint2D(x_max,y_current),aperture,cpt_fract);
		AddFracture(frac_mesh,net_mesh);
	}
	// addition of the vertical fractures
	for (int i=0;i<nb_division-1;i++){
		x_current+=fract_spacing_x;
		cpt_fract++;
		frac_mesh=FractureMesh(FluxPoint2D(x_current,y_min),FluxPoint2D(x_current,y_max),aperture,cpt_fract);
		AddFracture(frac_mesh,net_mesh);
	}
	// update the number of fractures
	net_mesh.cpt_fract=cpt_fract;
}

// Divide randomly few blocks of the domain and return the list of the divided blocks
vector<Domain> DivideSierpinskiDomain(pointcpp<double> domain_min,int nb_division,double fract_spacing_x,double fract_spacing_y,double aperture,NetworkMeshes& net_mesh,double Lx){
	// Parameters
	set<int> index_domain_list;int index;pair<int,int> indices;
	vector<Domain> Domain_list;
	// for each domain to divide
	for (int j=0;j<nb_division;j++){
		// Pick the domain to divide
		index=(int)(rand()/(double)RAND_MAX*nb_division*nb_division);
		while (index_domain_list.find(index)!=index_domain_list.end()){
			index=(int)(rand()/(double)RAND_MAX*nb_division*nb_division);
		}
		index_domain_list.insert(index);
		// determine its extremities
		indices=return_indices(index,nb_division);
		pointcpp<double> pt_min_current(domain_min.i+indices.first*fract_spacing_x,domain_min.j+indices.second*fract_spacing_y),
				pt_max_current(domain_min.i+(indices.first+1)*fract_spacing_x,domain_min.j+(indices.second+1)*fract_spacing_y);
		Domain_list.push_back(Domain(pt_min_current,pt_max_current));
		// divide it
		DivideDomain(pt_min_current,pt_max_current,nb_division,aperture,net_mesh,Lx);
	}
	return Domain_list;
}

void ComputeIntersections(NetworkMeshes& net_mesh){
	// 1. Variables
	bool intersection=true;
	FractureMesh fract1,fract2,new_fract1,new_fract2,new_fract3,new_fract4;
	Segment2D seg1,seg2;CgalPoint2D inter;FluxPoint2D pt_inter;
	int cpt_inter=net_mesh.cpt_inter,cpt_current;
	// 2. Loop on fracture mesh until there is no undefined intersection
	while (intersection){
		intersection=false;
		for (size_t i=0;i<net_mesh.meshes.size();i++){
			// 2.1. First segment definition
			fract1=net_mesh.meshes[i];
			seg1=fract1.ReturnSegment();
			for (size_t j=0;j<net_mesh.meshes.size();j++){
				bool new_node=false;
				if (i!=j){
					// 2.2. Second segment definition
					fract2=net_mesh.meshes[j];
					seg2=fract2.ReturnSegment();
					// 2.3. Check if there is intersection
					CGAL::Object result=CGAL::intersection(seg1,seg2);
					// if there is intersection
					if (!seg1.IdenticExtremities(seg2.source())&&!seg1.IdenticExtremities(seg2.target())&&CGAL::assign(inter,result)){
						// intersection definition
						cpt_current=net_mesh.ReturnIndex(inter);
						if (cpt_current==-1){cpt_current=cpt_inter+1;new_node=true;}
						pt_inter=FluxPoint2D(inter,cpt_current);
						// define the potential new subsegments from the first and second fracture
						new_fract1=FractureMesh(fract1.p_ori,pt_inter,fract1.aperture,fract1.fracture_index);
						new_fract2=FractureMesh(pt_inter,fract1.p_tar,fract1.aperture,fract1.fracture_index);
						new_fract3=FractureMesh(fract2.p_ori,pt_inter,fract2.aperture,fract2.fracture_index);
						new_fract4=FractureMesh(pt_inter,fract2.p_tar,fract2.aperture,fract2.fracture_index);
						// check if it is a real new intersection
						int cpt_true_inter=0;
						if (new_fract1.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (new_fract2.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (new_fract3.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (new_fract4.ReturnLength()>EPSILON){cpt_true_inter++;}
						if (cpt_true_inter>=3){intersection=true;}
						// if it is -> addition of the new subsegments
						if (intersection){
							if (new_fract1.ReturnLength()>EPSILON){
								net_mesh.meshes[i]=new_fract1;
								if (new_fract2.ReturnLength()>EPSILON){net_mesh.meshes.push_back(new_fract2);}
							}
							else{net_mesh.meshes[i]=new_fract2;}
							if (new_fract3.ReturnLength()>EPSILON){
								net_mesh.meshes[j]=new_fract3;
								if (new_fract4.ReturnLength()>EPSILON){net_mesh.meshes.push_back(new_fract4);}
							}
							else{net_mesh.meshes[j]=new_fract4;}
							if (new_node){cpt_inter++;}
							//net_mesh.print_DFN();
							break;
						}
					}
				}
			}
			if (intersection){break;}
		}
	}
	// 3. Update network values
	net_mesh.cpt_inter=cpt_inter;
}


// Define Sierpinski structures (whole horizontal and vertical fractures)
void NetworkMeshesGenerationSierpinski(string code_path, string file_name_domain,NetworkMeshes& net_mesh){
	// 0. Variables
	pointcpp<double> pt_min_init=net_mesh.domain.min_pt,pt_max_init=net_mesh.domain.max_pt;
	if (net_mesh.domain.domain_size_x()!=net_mesh.domain.domain_size_y()){cout << "WARNING in NetworkMeshesGenerationSierpinski (NetworkMeshes.cpp): option not implemented" << endl;}
	double Lx=net_mesh.domain.domain_size_x();
	// parameters
	int nb_division,level_division,seed; double aperture; string generation_option;
	string file_name = code_path + "/Input/DFN_files/";
	file_name += file_name_domain;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> generation_option >> nb_division >> level_division >> aperture >> seed;
	}
	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;
	}
	fichier.close();
	// 1. First domain division
	srand(seed);
	DivideDomain(pt_min_init,pt_max_init,nb_division,aperture,net_mesh,Lx);
	// extra-fracture on the top of the domain
	/*FractureMesh frac_mesh(FluxPoint2D(pt_min_init.i,pt_max_init.j),FluxPoint2D(pt_max_init.i,pt_max_init.j),aperture);
	AddFracture(frac_mesh,net_mesh);*/
	// 2. Next domain division
	vector<Domain> divided_domains;divided_domains.push_back(net_mesh.domain);
	double fract_spacing_x=(net_mesh.domain.max_pt.i-net_mesh.domain.min_pt.i)/nb_division,fract_spacing_y=(net_mesh.domain.max_pt.j-net_mesh.domain.min_pt.j)/nb_division;
	net_mesh.max_fract_spacing=fract_spacing_y;
	for (int i=1;i<level_division;i++){
		// division of the domain
		vector<Domain> new_divided_domains,current_divided_domains;
		for (vector<Domain>::iterator it1=divided_domains.begin();it1!=divided_domains.end();it1++){
			fract_spacing_x=(it1->max_pt.i-it1->min_pt.i)/nb_division,fract_spacing_y=(it1->max_pt.j-it1->min_pt.j)/nb_division;
			current_divided_domains=DivideSierpinskiDomain(it1->min_pt,nb_division,fract_spacing_x,fract_spacing_y,aperture,net_mesh,Lx);
			new_divided_domains.insert(new_divided_domains.end(),current_divided_domains.begin(),current_divided_domains.end());
		}
		divided_domains=new_divided_domains;
	}
	// 3. Domain translation to work on a periodic domain in the vertical direction
	//net_mesh=Vertical_Translation(net_mesh,net_mesh.max_fract_spacing*0.5);
	// 4. Intersections determination
	ComputeIntersections(net_mesh);
}

/******************* FUNCTION TO GENERATE A REALISTIC FRACTURE NETWORK **********************/
void NetworkMeshesGenerationRealistic(string code_path, string file_name_domain,NetworkMeshes& net_mesh){
	// 0. Variables
	if (net_mesh.domain.domain_size_x()!=net_mesh.domain.domain_size_y()){cout << "WARNING in NetworkMeshesGenerationRealistic (NetworkMeshes.cpp): option not implemented" << endl;}
	// parameters
	double C_density,D_dim,coeff_theta,fract_aperture,r_min;int seed=0;
	string file_name = code_path + "/Input/DFN_files/",generation_option;
	file_name += file_name_domain;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> generation_option >> C_density >> D_dim >> fract_aperture >> r_min >> coeff_theta >> seed;
	}
	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;
	}
	fichier.close();

	// nb_fract is equal to the integer of i_fract rounded up.
	int nb_fract=std::ceil(C_density*std::pow(r_min,-(double)D_dim));
	int cpt_fract=0;

	FractureMesh frac_mesh; Segment2D seg;double fract_length,angle;CgalPoint2D fract_center;
	double Lx=net_mesh.domain.domain_size_x(),Ly=net_mesh.domain.domain_size_y();

	srand(seed);
	while (cpt_fract<nb_fract){
		cpt_fract++;
		// center of the fracture
		fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
		// length of the fracture
		fract_length=std::pow(C_density/(double)cpt_fract,1./(double)D_dim)*Lx;
		//angle of the fracture
		if ((double)cpt_fract/2.==(int)cpt_fract/2){angle=coeff_theta*PI/180;}
		else {angle=(coeff_theta+120)*PI/180;}
		// definition of a segment
		seg=Segment2D(fract_center,fract_length,angle);
		net_mesh.domain.SegmentIntersectDomain(seg,seg);
		frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
		AddFracture(frac_mesh,net_mesh);
	}
	// 4. Intersections determination
	ComputeIntersections(net_mesh);
}

void NetworkMeshesGenerationRealistic2(string code_path, string file_name_domain,NetworkMeshes& net_mesh){
	// 0. Variables
	if (net_mesh.domain.domain_size_x()!=net_mesh.domain.domain_size_y()){cout << "WARNING in NetworkMeshesGenerationRealistic (NetworkMeshes.cpp): option not implemented" << endl;}
	// parameters
	double C_density,D_dim,coeff_theta1,coeff_theta2,fract_aperture,r_min;int seed=0;
	string file_name = code_path + "/Input/DFN_files/",generation_option;
	file_name += file_name_domain;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		fichier >> generation_option >> C_density >> D_dim >> fract_aperture >> r_min >> coeff_theta1 >> coeff_theta2 >> seed;
	}
	else{
		cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;
	}
	fichier.close();

	// nb_fract is equal to the integer of i_fract rounded up.
	int nb_fract=std::ceil(C_density*std::pow(r_min,-(double)D_dim));
	int cpt_fract=0;

	FractureMesh frac_mesh; Segment2D seg;double fract_length,angle;CgalPoint2D fract_center;
	double Lx=net_mesh.domain.domain_size_x(),Ly=net_mesh.domain.domain_size_y();

	srand(seed);
	while (cpt_fract<nb_fract){
		cpt_fract++;
		// center of the fracture
		fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
		// length of the fracture
		fract_length=std::pow(C_density/(double)cpt_fract,1./(double)D_dim)*Lx;
		//angle of the fracture
		if ((double)cpt_fract/2.==(int)cpt_fract/2){angle=coeff_theta1*PI/180;}
		else {angle=coeff_theta2*PI/180;}
		// definition of a segment
		seg=Segment2D(fract_center,fract_length,angle);
		net_mesh.domain.SegmentIntersectDomain(seg,seg);
		frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
		AddFracture(frac_mesh,net_mesh);
	}
	// 4. Intersections determination
	ComputeIntersections(net_mesh);
}

void NetworkMeshesGenerationRealistic3(string code_path, string file_name_domain,NetworkMeshes& net_mesh){
        // 0. Variables
        if (net_mesh.domain.domain_size_x()!=net_mesh.domain.domain_size_y()){cout << "WARNING in NetworkMeshesGenerationRealistic (NetworkMeshes.cpp): option not implemented" << endl;}
        // parameters
        double C_density,D_dim,fract_aperture,r_min;int seed=0;
	double b_min, b_max, mean_lnb, RSD_lnb;
        string file_name = code_path + "/Input/DFN_files/",generation_option;
        file_name += file_name_domain;
        ifstream fichier(file_name.c_str());
        if (fichier.is_open()){
                //fichier >> generation_option >> C_density >> D_dim >> fract_aperture >> r_min >> seed;
        	fichier >> generation_option >> C_density >> D_dim >> b_min >> b_max >> mean_lnb >> RSD_lnb >> r_min >> seed;
	}
        else{
                cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;
        }
        fichier.close();

        // nb_fract is equal to the integer of i_fract rounded up.
        int nb_fract=std::ceil(C_density*std::pow(r_min,-(double)D_dim));
        int cpt_fract=0;

        FractureMesh frac_mesh; Segment2D seg;double fract_length,angle;CgalPoint2D fract_center;
        double Lx=net_mesh.domain.domain_size_x(),Ly=net_mesh.domain.domain_size_y();

        srand(seed);
        while (cpt_fract<nb_fract){
                cpt_fract++;
                // center of the fracture
                fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
                // length of the fracture
                fract_length=std::pow(C_density/(double)cpt_fract,1./(double)D_dim)*Lx;
                //angle of the fracture
                angle=Uniform()*PI;
cout << "angle = " << angle << endl;
                // definition of a segment
                seg=Segment2D(fract_center,fract_length,angle);
                net_mesh.domain.SegmentIntersectDomain(seg,seg);
		//fract_aperture=b_min;
		fract_aperture=exp(sqrt(2)*RSD_lnb*boost::math::erf_inv(Uniform()*(erf((log(b_max)-mean_lnb)/(sqrt(2)*RSD_lnb))-erf((log(b_min)-mean_lnb)/(sqrt(2)*RSD_lnb)))+erf((log(b_min)-mean_lnb)/(sqrt(2)*RSD_lnb)))+mean_lnb);
		//fract_aperture=boost::math::erf_inv(b_max);
                frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
                AddFracture(frac_mesh,net_mesh);
        }
        // 4. Intersections determination
        ComputeIntersections(net_mesh);
}

NetworkMeshes::NetworkMeshes(double density_param, double exponent_param, Parameters param, Domain domain_){
	// 1. Initial properties
	domain=domain_;
        cpt_inter=-1;cpt_mesh=-1;cpt_fract=-1;
	// 2. Network generation
	int cpt_fract=0;double density=0.0,rnd_value;
	double fract_aperture;
	FractureMesh frac_mesh; Segment2D seg;double fract_length,angle;CgalPoint2D fract_center;
	double Lx=domain.domain_size_x(),Ly=domain.domain_size_y();
	srand(param.seed_simu);
	if (param.generation_option_DFN=="generation_realistic2"){	// Zitong's study: C,D for the length, two angles, constant aperture
		// nb_fract is equal to the integer of i_fract rounded up.
		int nb_fract=std::ceil(density_param*std::pow(param.r_min,-(double)exponent_param));
		while (cpt_fract<nb_fract){
			cpt_fract++;
			fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);// center of the fracture
		 	fract_length=std::pow(density_param/(double)cpt_fract,1./(double)exponent_param)*Lx;// length of the fracture
			//angle of the fracture
			if ((double)cpt_fract/2.==(int)cpt_fract/2){angle=param.coeff_theta1*PI/180;}
			else {angle=param.coeff_theta2*PI/180;}
			// definition of a segment
			seg=Segment2D(fract_center,fract_length,angle);
			domain.SegmentIntersectDomain(seg,seg);
			frac_mesh=FractureMesh(seg,param.fract_aperture,cpt_fract-1);
			AddFracture(frac_mesh,*this);
		}
	}
	else if (param.generation_option_DFN=="generation_realistic3"){		// Guofeng's study: a,p for the length, random angles, distributed aperture, correlated length
		while (density<density_param){
                        fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);// center of the fracture
                       double rnd_value=Uniform();
                       fract_length=pow(rnd_value,1./(1-exponent_param));
			//angle of the fracture
                	angle=Uniform()*PI;
                        // definition of a segment
                        seg=Segment2D(fract_center,fract_length,angle);
                        domain.SegmentIntersectDomain(seg,seg);
			rnd_value=1-rnd_value;
			fract_aperture=exp(sqrt(2)*param.RSD_lnb*boost::math::erf_inv(rnd_value*(erf((log(param.b_max)-param.mean_lnb)/(sqrt(2)*param.RSD_lnb))-erf((log(param.b_min)-param.mean_lnb)/(sqrt(2)*param.RSD_lnb)))+erf((log(param.b_min)-param.mean_lnb)/(sqrt(2)*param.RSD_lnb)))+param.mean_lnb);
                        frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
                        AddFracture(frac_mesh,*this);
			density+=seg.length*seg.length/(Lx*Ly);
			cout << "length = " << fract_length << ", aperture = " << fract_aperture << endl;
        	}
	}
	else if (param.generation_option_DFN=="generation_realistic4"){		// Guofeng's study with not correlated length and aperture
                while (density<density_param){
                        fract_center=CgalPoint2D(-0.5*Lx+Uniform()*Lx,-0.5*Ly+Uniform()*Ly);
			fract_length=pow(Uniform(),1./(1-exponent_param));
                        angle=Uniform()*PI;
			seg=Segment2D(fract_center,fract_length,angle);
                        domain.SegmentIntersectDomain(seg,seg);
			fract_aperture=exp(sqrt(2)*param.RSD_lnb*boost::math::erf_inv(Uniform()*(erf((log(param.b_max)-param.mean_lnb)/(sqrt(2)*param.RSD_lnb))-erf((log(param.b_min)-param.mean_lnb)/(sqrt(2)*param.RSD_lnb)))+erf((log(param.b_min)-param.mean_lnb)/(sqrt(2)*param.RSD_lnb)))+param.mean_lnb);                      
                        frac_mesh=FractureMesh(seg,fract_aperture,cpt_fract-1);
                        AddFracture(frac_mesh,*this);
			density+=seg.length*seg.length/(Lx*Ly);
                }
        }
	else{cout << "WARNING in NetworkMeshes.cpp: generation option not implemented" << endl;}
	// Intersections determination
	ComputeIntersections(*this);
}

/****************** FUNCTION TO READ FRACTURE NETWORK DESCRIPTION **************************/
void NetworkMeshesFromFile(string code_path, string file_name_domain,NetworkMeshes& net_mesh){
	// local variables
	string str_line,x_str,y_str;
	int pos1,pos2,nb_segment,cpt_fract=0;
	double aperture;string generation_option;
	// Domain definition from files
	string file_name = code_path + "/Input/DFN_files/";
	file_name += file_name_domain;FractureMesh frac_mesh;Segment2D seg;
	ifstream fichier(file_name.c_str());
	if (fichier.is_open()){
		// number of segments to describe the fracture network
		fichier >> generation_option >> nb_segment;
		net_mesh.meshes.reserve(nb_segment);
		// each line describes a segment (or mesh) composing the fracture network
		double x1=0.0,y1=0.0,x2=0.0,y2=0.0;
		for (int i=0;i<nb_segment;i++){
			if (!(fichier >> x1 >> y1 >> x2 >> y2 >> aperture)){break;}
			CgalPoint2D p_origin(x1,y1);
			CgalPoint2D p_target(x2,y2);
			seg=Segment2D(p_origin,p_target);
			frac_mesh=FractureMesh(seg,aperture,cpt_fract);
			cpt_fract++;
			AddFracture(frac_mesh,net_mesh);
		}
	}
	else{cout << "WARNING in Parameters (Parameters.cpp) : unknown file" << endl;}
	fichier.close();
	// 4. Intersections determination
	ComputeIntersections(net_mesh);
}

FractureMesh NetworkMeshes::return_mesh(int mesh_ind){
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->mesh_index == mesh_ind){
			return *it;
		}
	}
	cout << "WARNING in return_mesh (NetworkMeshes.cpp): mesh number " << mesh_ind << "not found" << endl;
	return FractureMesh();
}

FractureMesh NetworkMeshes::return_mesh(int mesh_ind) const{
	for (vector<FractureMesh>::const_iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->mesh_index == mesh_ind){
			return *it;
		}
	}
	cout << "WARNING in return_mesh (NetworkMeshes.cpp): mesh number " << mesh_ind << "not found" << endl;
	return FractureMesh();
}

void NetworkMeshes::print_DFN(){
	cout << endl;
	cout << "DFN with " << meshes.size() << " meshes:" << endl;
	for (int i=0;i<meshes.size();i++){
		cout << "(" << meshes[i].p_ori.p.x() << "," << meshes[i].p_ori.p.y() << ") " << meshes[i].p_ori.index << " (" << meshes[i].p_tar.p.x() << "," << meshes[i].p_tar.p.y() << ") "
		<< meshes[i].p_tar.index << " " << meshes[i].velocity << " " << meshes[i].aperture << " " << meshes[i].mesh_index << " ";
		for (int j=0;j<meshes[i].neigh_meshes.size();j++){
			cout << meshes[i].neigh_meshes[j] << " ";
		}
		cout << meshes[i].fracture_index;
		cout << endl;
	}
	cout << endl;
}

void NetworkMeshes::print_fractures(){
	cout << endl;
	cout << fractures.size() << " fractures:" << endl;
	for (FracturesMap::iterator it=fractures.begin();it!=fractures.end();it++){
		cout << "fracture index = " << it->first << ", origin index = " << it->second.first.index << " (" << it->second.first.p.x() << "," << it->second.first.p.y() << ")" << ", target index = " << it->second.second.index
				<< " (" << it->second.second.p.x() << "," << it->second.second.p.y() << ")"  << endl;
	}
	cout << endl;
}

// This prints the mesh coordinates Feb 2018
void NetworkMeshes::print_DFN_in_file(Parameters param){
	string file_name=param.code_path+"/Output/DFN.txt";
	ofstream output(file_name.c_str(),ofstream::out);
	for (int i=0;i<meshes.size();i++){
		output << meshes[i].p_ori.p.x() << " " << meshes[i].p_ori.p.y() << " " << meshes[i].p_tar.p.x() << " " << meshes[i].p_tar.p.y() << " " << meshes[i].aperture << " " << meshes[i].velocity << endl;
	}
	output.close();
}

void NetworkMeshes::print_DFN_in_file(string output_path,int i){
	stringstream ss; ss << i;
	string file_name=output_path+"/dfn/DFN"+ss.str()+".txt";
	ofstream output(file_name.c_str(),ofstream::out);
	for (int i=0;i<meshes.size();i++){
		output << meshes[i].p_ori.p.x() << " " << meshes[i].p_ori.p.y() << " " << meshes[i].p_tar.p.x() << " " << meshes[i].p_tar.p.y() << endl;
	}
	output.close();
}

void NetworkMeshes::print_DFN_in_file_initial(string output_path,int i){
        stringstream ss; ss << i;
        string file_name=output_path+"/dfn_initial/DFN"+ss.str()+".txt";
        ofstream output(file_name.c_str(),ofstream::out);
        for (int i=0;i<meshes.size();i++){
                output << meshes[i].p_ori.p.x() << " " << meshes[i].p_ori.p.y() << " " << meshes[i].p_tar.p.x() << " " << meshes[i].p_tar.p.y() << " " << meshes[i].aperture << endl;
        }
        output.close();
}

void NetworkMeshes::print_DFN_in_file_initial(Parameters param){
        string file_name=param.code_path+"/Output/DFN_init.txt";
	ofstream output(file_name.c_str(),ofstream::out);
	for (int i=0;i<meshes.size();i++){
		output << meshes[i].p_ori.p.x() << " " << meshes[i].p_ori.p.y() << " " << meshes[i].p_tar.p.x() << " " << meshes[i].p_tar.p.y() << " " << meshes[i].aperture << " " << meshes[i].velocity << endl;
	}
	output.close();
}

void NetworkMeshes::print_DFN_in_file_with_aperture_delta(Parameters param,const NetworkMeshes& initial_mesh){
	string file_name=param.code_path+"/Output/DFN_aperture_delta.txt";
	ofstream output(file_name.c_str(),ofstream::out);
	std::map<int,double> init_aperture;
	for (size_t i=0;i<initial_mesh.meshes.size();i++){
		init_aperture[initial_mesh.meshes[i].mesh_index]=initial_mesh.meshes[i].aperture;
	}
	for (int i=0;i<meshes.size();i++){
		double b_init = 0.0;
		std::map<int,double>::iterator it = init_aperture.find(meshes[i].mesh_index);
		if (it!=init_aperture.end()){b_init = it->second;}
		double delta_b = meshes[i].aperture - b_init;
		output << meshes[i].p_ori.p.x() << " " << meshes[i].p_ori.p.y() << " "
		       << meshes[i].p_tar.p.x() << " " << meshes[i].p_tar.p.y() << " "
		       << meshes[i].aperture << " " << meshes[i].velocity << " " << delta_b << endl;
	}
	output.close();
}

// print aperture data 
//void NetworkMeshes::print_DFN_in_file_aperture(string output_path,int i){
  //       stringstream ss; ss << i;
    //     string file_name=output_path+"/dfn_aperture/DFN"+ss.str()+".txt";
      //   ofstream output(file_name.c_str(),ofstream::out);
       //  for (int i=0;i<meshes.size();i++){
        //       	output << meshes[i].aperture << endl;
	  //}
        // output.close();
//}


// return the map of node index and a list of their corresponding fractures
NodeFracturesMap NetworkMeshes::ReturnNodeFracturesMap(){
	// 0. Variables
	NodeFracturesMap node_fract_map;
	// 1. Loop on the fractures (or meshes)
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		// 1.1. Study of the first extremity
		NodeInsertionMap(it->p_ori.index,it->p_tar.index,*it,node_fract_map);
		// 1.2. Study of the second extremity  (if not defined already)
		NodeInsertionMap(it->p_tar.index,it->p_ori.index,*it,node_fract_map);
	}
	return node_fract_map;
}

// return the list of nodes connected to the domain borders
set<int> NetworkMeshes::return_connected_nodes(){
	set<int> connected_nodes;NodeFracturesMap node_fract_map=ReturnNodeFracturesMap();
	for (BordersMap::iterator it=border_map.begin();it!=border_map.end();it++){
		connected_nodes.insert(it->first);
		add_neigh_nodes(connected_nodes,node_fract_map,it->first);
	}
	return connected_nodes;
}

vector<FractureMesh> NetworkMeshes::return_mesh_from_node_ori(int node_index){
	vector<FractureMesh> mesh_list;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (it->p_ori.index==node_index){
			mesh_list.push_back(*it);
		}
	}
	return mesh_list;
}

// considering that the index of each mesh and of the neighbors of each mesh are redefined after this function
NetworkMeshes NetworkMeshes::remove_fractures_from_velocity(Parameters param,double velocity_threshold,string option){
	NetworkMeshes new_net_mesh=*this;
	std::vector<FractureMesh> new_meshes;
	ComputeFlowVelocities(*this,param,domain,option);
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
			if (fabs(it->velocity)>velocity_threshold){
				new_meshes.push_back(*it);
			}
	}
	new_net_mesh.meshes=new_meshes;
	return new_net_mesh;
}

// return a new networkmeshes from a list of nodes
NetworkMeshes NetworkMeshes::return_network_from_nodes(set<int> connected_nodes){
	NetworkMeshes new_net_mesh;
	// Domain
	new_net_mesh.domain=domain;new_net_mesh.cpt_fract=cpt_fract;
	// Connected node numbering
	map<int,int> node_numbering;//<new_index,old_index>
	int cpt_node=-1;

	for (set<int>::iterator it1=connected_nodes.begin();it1!=connected_nodes.end();it1++){
		cpt_node++;
		node_numbering[cpt_node]=*it1;
	}
	// nodes and mesh selection
	vector<FractureMesh> mesh_list;FractureMesh mesh;int old_index,new_index;
	int cpt_mesh=0;map<int,int>::iterator it_map;
	for (map<int,int>::iterator it1=node_numbering.begin();it1!=node_numbering.end();it1++){
		new_index=it1->first;old_index=it1->second;
		// Border map
		if (border_map.find(old_index)!=border_map.end()){new_net_mesh.border_map[new_index]=border_map[old_index];}
		// Node map
		new_net_mesh.nodes_map[new_index]=nodes_map[old_index];
		// Inter list
		if (inter_list.find(old_index)!=inter_list.end()){new_net_mesh.inter_list.insert(new_index);}
		// Meshes
		mesh_list=return_mesh_from_node_ori(old_index);
		for (vector<FractureMesh>::iterator it2=mesh_list.begin();it2!=mesh_list.end();it2++){
			if (connected_nodes.find(it2->p_ori.index)!=connected_nodes.end()&&connected_nodes.find(it2->p_tar.index)!=connected_nodes.end()){
				mesh=*it2;
				for (map<int,int>::iterator it3=node_numbering.begin();it3!=node_numbering.end();it3++){
					if (it3->second==mesh.p_ori.index){mesh.p_ori.index=it3->first;}
				}
				for (map<int,int>::iterator it3=node_numbering.begin();it3!=node_numbering.end();it3++){
					if (it3->second==mesh.p_tar.index){mesh.p_tar.index=it3->first;}
				}
				mesh.mesh_index=cpt_mesh;
				new_net_mesh.meshes.push_back(mesh);
				cpt_mesh++;
			}
		}
	}

	new_net_mesh.cpt_inter=cpt_node;
	return new_net_mesh;
}

// return the nodes corresponding to meshes with a flow velocity larger than a given threshold
set<int> NetworkMeshes::return_nodes_from_velocity(Parameters param,double velocity_threshold,string option){
	set<int> list_nodes;
	// compute flow velocities
	ComputeFlowVelocities(*this,param,domain,option);
	// selection of the meshes with a flow velocity larger than the velocity threshold
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (fabs(it->velocity)>velocity_threshold){
			list_nodes.insert(it->p_ori.index);
			list_nodes.insert(it->p_tar.index);
		}
	}
	return list_nodes;
}

NetworkMeshes NetworkMeshes::return_backbone(Parameters param,double velocity_threshold,string option){
	NetworkMeshes new_network_meshes;
	// selection of the network with nodes connected to the borders
	set<int> connected_nodes=return_connected_nodes();
	if (connected_nodes.size()==0){return new_network_meshes;}
	new_network_meshes=return_network_from_nodes(connected_nodes);
	// selection of the meshes with flow velocities larger than a given threshold
	set<int> nodes_from_velocity=new_network_meshes.return_nodes_from_velocity(param,velocity_threshold,option);
	new_network_meshes=new_network_meshes.return_network_from_nodes(nodes_from_velocity);
	if (new_network_meshes.meshes.size()==0){return new_network_meshes;}
	new_network_meshes=new_network_meshes.remove_fractures_from_velocity(param,velocity_threshold,option);
	return new_network_meshes;
}

void NetworkMeshes::EvaluateFlowVelocities(ublas_vector DFN_potential){
	FractureMesh fract_mesh;double velocity;
	int nb_nodes = DFN_potential.size();
	for (int i=0;i<meshes.size();i++){
		fract_mesh=meshes[i];
		if (meshes[i].aperture<=0.0){
			meshes[i].velocity=0.0;
			continue;
		}
		if (fract_mesh.p_ori.index<0 || fract_mesh.p_ori.index>=nb_nodes ||
		    fract_mesh.p_tar.index<0 || fract_mesh.p_tar.index>=nb_nodes){
			cout << "WARNING in EvaluateFlowVelocities (NetworkMeshes.cpp): node index out of range"
			     << " mesh_index=" << meshes[i].mesh_index
			     << " p_ori=" << fract_mesh.p_ori.index
			     << " p_tar=" << fract_mesh.p_tar.index
			     << " nb_nodes=" << nb_nodes << endl;
			meshes[i].velocity=0.0;
			continue;
		}
		double length = fract_mesh.ReturnLength();
		double head_ori = DFN_potential(fract_mesh.p_ori.index);
		double head_tar = DFN_potential(fract_mesh.p_tar.index);
		velocity=meshes[i].ReturnConductance()*(head_ori-head_tar)/length;
		if (!std::isfinite(velocity) || !std::isfinite(head_ori) || !std::isfinite(head_tar) ||
		    !std::isfinite(length) || length<=0.0){
			cout << "WARNING in EvaluateFlowVelocities (NetworkMeshes.cpp): invalid flow values"
			     << " mesh_index=" << meshes[i].mesh_index
			     << " p_ori=" << fract_mesh.p_ori.index
			     << " p_tar=" << fract_mesh.p_tar.index
			     << " head_ori=" << head_ori
			     << " head_tar=" << head_tar
			     << " length=" << length
			     << " aperture=" << meshes[i].aperture
			     << " velocity=" << velocity << endl;
		}
		meshes[i].velocity=velocity;
	}
}

pair<double,double> NetworkMeshes::return_min_max_vel(){
	if (meshes.empty()){
		cout << "WARNING in return_min_max_vel (NetworkMeshes.cpp): meshes is empty" << endl;
		return make_pair(-1,-1);
	}
	pair<double,double> min_max_vel=make_pair(abs(meshes.begin()->velocity),abs(meshes.begin()->velocity));
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		if (abs(it->velocity)<min_max_vel.first){min_max_vel.first=abs(it->velocity);}
		if (abs(it->velocity)>min_max_vel.second){min_max_vel.second=abs(it->velocity);}
	}
	return min_max_vel;
}

double NetworkMeshes::return_ave_vel(){
	if (meshes.empty()){
		cout << "WARNING in return_ave_vel (NetworkMeshes.cpp): meshes is empty" << endl;
		return -1;
	}
	double ave_vel=0;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		ave_vel+=abs(it->velocity);
	}
	ave_vel = ave_vel/meshes.size();
	return ave_vel;
}

// Define the map fractures: index and extremities of each fracture with positive flow velocity going from the origin to the end of the fractures
bool NetworkMeshes::define_fracture_map(){
	// Variables
	FluxPoint2D ext1,ext2; bool p_ori_def_twice,p_tar_def_twice,ext1_define,ext2_define,test=true;

	cout << "Define extremities" << endl;
	int new_cpt_fract=-1;int current_node_index;
	map<int, set<int> > fract_mesh_indices;	// map<fracture index,set<mesh indices> >
	// Define the extremities of each fracture (not defined yet, as some pieces disappear when the backbone is computed)
	for (set<int>::iterator it_fract=fract_indices.begin();it_fract!=fract_indices.end();it_fract++){
		set<int> nodes_defined_once_index,nodes_defined_once_treated_index;
		ext1_define=false;ext2_define=false;
		for (vector<FractureMesh>::iterator it1=meshes.begin();it1!=meshes.end();it1++){
			p_ori_def_twice=false;p_tar_def_twice=false;
			if (it1->fracture_index==*it_fract){
				for (vector<FractureMesh>::iterator it2=meshes.begin();it2!=meshes.end();it2++){
					if (it1!=it2&&it2->fracture_index==*it_fract){
						if ((it1->p_ori.index==it2->p_ori.index)||(it1->p_ori.index==it2->p_tar.index)){p_ori_def_twice=true;}
						if ((it1->p_tar.index==it2->p_ori.index)||(it1->p_tar.index==it2->p_tar.index)){p_tar_def_twice=true;}
					}
				}
				if (!p_ori_def_twice){nodes_defined_once_index.insert(it1->p_ori.index);}
				if (!p_tar_def_twice){nodes_defined_once_index.insert(it1->p_tar.index);}
			}
		}
		// Define the corresponding fractures (several fractures for discontinuous cases)
		for (set<int>::iterator it_node=nodes_defined_once_index.begin();it_node!=nodes_defined_once_index.end();it_node++){
			if (nodes_defined_once_treated_index.find(*it_node)==nodes_defined_once_treated_index.end()){		// node not treated yet
				set<int> mesh_indices;new_cpt_fract++;
				current_node_index=*it_node;bool final_node_not_found=true;
				while (final_node_not_found){
					for (vector<FractureMesh>::iterator it1=meshes.begin();it1!=meshes.end();it1++){
						if (it1->fracture_index==*it_fract&&(it1->p_ori.index==current_node_index||it1->p_tar.index==current_node_index)){
							mesh_indices.insert(it1->mesh_index);
							if (it1->p_ori.index==current_node_index){current_node_index=it1->p_tar.index;}
							else if (it1->p_tar.index==current_node_index){current_node_index=it1->p_ori.index;}
							if(nodes_defined_once_index.find(current_node_index)!=nodes_defined_once_index.end()){		// node is found
								final_node_not_found=false;break;
							}
						}
					}
				}
				bool ext1_define=false,ext2_define=false;
				for (vector<FractureMesh>::iterator it_mesh=meshes.begin();it_mesh!=meshes.end();it_mesh++){
					if (it_mesh->p_ori.index==*it_node){ext1=it_mesh->p_ori;ext1_define=true;}
					else if (it_mesh->p_tar.index==*it_node){ext1=it_mesh->p_tar;ext1_define=true;}
					if (it_mesh->p_ori.index==current_node_index){ext2=it_mesh->p_ori;ext2_define=true;}
					else if (it_mesh->p_tar.index==current_node_index){ext2=it_mesh->p_tar;ext2_define=true;}
					if (ext1_define&&ext2_define){break;}
				}
				fractures[new_cpt_fract]=make_pair(ext1,ext2);
				fract_mesh_indices[new_cpt_fract]=mesh_indices;
				nodes_defined_once_treated_index.insert(*it_node);
				nodes_defined_once_treated_index.insert(current_node_index);
			}
		}
	}

	// Update fracture index in meshes
	for (map<int, set<int> >::iterator it1=fract_mesh_indices.begin();it1!=fract_mesh_indices.end();it1++){
		for (set<int>::iterator it2=it1->second.begin();it2!=it1->second.end();it2++){
			for (int i_mesh=0;i_mesh<meshes.size();i_mesh++){
				if (meshes[i_mesh].mesh_index==*it2){meshes[i_mesh].fracture_index=it1->first;}
			}
		}
	}

	//print_fractures();
	//print_DFN();

	cout << "Redefine fractures depending on their flow velocity" << endl;
	// Redefine fractures such as the flow velocity in each fracture goes from its origin to its target (will require to divide some fractures into 2 or more fractures)
	double previous_velocity=1,current_velocity=1;int cpt_fract_new=0;
	FluxPoint2D start_fract;
	FracturesMap new_fractures; map<int,int> mesh_fract_index;	// <mesh index,new fracture index>
	// For each fracture
	for (FracturesMap::iterator it_fract=fractures.begin();it_fract!=fractures.end();it_fract++){
		start_fract=it_fract->second.first;ext1=start_fract;ext2=start_fract;previous_velocity=1;current_velocity=1;
		// ext1 and ext2: extremities of the previously studied mesh
		// Loop on the corresponding meshes until reaching the second extremity of the fracture
		while (ext2.index!=it_fract->second.second.index){
			previous_velocity=current_velocity;

			// look for the mesh with the current node index and on the studied fracture
			for (vector<FractureMesh>::iterator it1=meshes.begin();it1!=meshes.end();it1++){
				if (it1->fracture_index==it_fract->first&&((it1->p_ori.index==ext2.index&&it1->p_tar.index!=ext1.index)||(it1->p_tar.index==ext2.index&&it1->p_ori.index!=ext1.index))){
					if (it1->p_ori.index==ext2.index){current_velocity=it1->velocity;ext1=ext2;ext2=it1->p_tar;}
					else if (it1->p_tar.index==ext2.index){current_velocity=-it1->velocity;ext1=ext2;ext2=it1->p_ori;}
					mesh_fract_index[it1->mesh_index]=cpt_fract_new;
					break;
				}
			}

			if (current_velocity*previous_velocity<0&&start_fract.index!=ext1.index){	// not at the first segment of the fracture
				if (previous_velocity<0){new_fractures[cpt_fract_new]=make_pair(ext1,start_fract);}		// we store the current fracture
				else{new_fractures[cpt_fract_new]=make_pair(start_fract,ext1);}
				cpt_fract_new++;
				start_fract=ext1;
				// look for the previous ext1 to not pick it again
				for (vector<FractureMesh>::iterator it1=meshes.begin();it1!=meshes.end();it1++){
					if (it1->fracture_index==it_fract->first&&((it1->p_ori.index==start_fract.index&&it1->p_tar.index!=ext2.index)||(it1->p_tar.index==start_fract.index&&it1->p_ori.index!=ext2.index))){
						if (it1->p_ori.index==start_fract.index){ext1=it1->p_tar;}
						else if (it1->p_tar.index==start_fract.index){ext1=it1->p_ori;}
					}
				}
				ext2=start_fract;previous_velocity=current_velocity;
			}

		}
		if (current_velocity<0){new_fractures[cpt_fract_new]=make_pair(ext2,start_fract);}		// we store the current fracture
		else{new_fractures[cpt_fract_new]=make_pair(start_fract,ext2);}
		cpt_fract_new++;
	}
	fractures=new_fractures;

	// Update fracture index in meshes
	for (map<int,int>::iterator it_new_index=mesh_fract_index.begin();it_new_index!=mesh_fract_index.end();it_new_index++){
		for (int i_mesh=0;i_mesh<meshes.size();i_mesh++){
			if (meshes[i_mesh].mesh_index==it_new_index->first){meshes[i_mesh].fracture_index=it_new_index->second;break;}
		}
	}

	//print_fractures();
	//print_DFN();

	cout << "Check positive flow velocity" << endl;
	// check that the flow velocity is positive from the origin to the end of the fractures for each fracture and all along the fractures
	int current_index;
	for (FracturesMap::iterator it_fract=fractures.begin();it_fract!=fractures.end();it_fract++){
		current_index=it_fract->second.first.index;
		while (current_index!=it_fract->second.second.index){
			for (vector<FractureMesh>::iterator it1=meshes.begin();it1!=meshes.end();it1++){
				if (current_index==it1->p_ori.index&&it1->fracture_index==it_fract->first&&it1->velocity>0){
					current_index=it1->p_tar.index;
					break;
				}
				else if (current_index==it1->p_tar.index&&it1->fracture_index==it_fract->first&&it1->velocity<0){
					current_index=it1->p_ori.index;
					break;
				}
			}
		}
	}
	cout << "end check" << endl;
	return test;
}


// Compute and return the volumetric fracture density
double NetworkMeshes::ReturnVolumetricFractureDensity(){
	double vol_dens=0.0;
	for (vector<FractureMesh>::iterator it=meshes.begin();it!=meshes.end();it++){
		//VRG:24 jan 2018 replacing vol_dens (which is teachnical a volumetric fracture frequency) with an actual fracture density as defined by Renshaw.
		//vol_dens+=  (it->ReturnLength()/2)*(it->ReturnLength()/2);
		vol_dens+=it->aperture*it->ReturnLength();
		//cout << "Return Length0" << it->ReturnLength() << endl;
	}
	return vol_dens/(domain.domain_size_x()*domain.domain_size_y());
}

/****** FUNCTIONS FOR THE BACKBONE DEFINITION ***************/
// addition of the node in the numbering map
void NodeInsertionMap(int node1,int node2,FractureMesh fract_mesh,NodeFracturesMap& node_fract_map){
	HydrauPropMap prop_map;
	// if already defined
	NodeFracturesMap::iterator it1=node_fract_map.find(node1);
	if (it1!=node_fract_map.end()){prop_map=it1->second;}
	// addition of the studied segment information in the submap
	prop_map[node2]=HydraulicProperties(fract_mesh);
	// addition of the information in the main map
	node_fract_map[node1]=prop_map;
}

void add_neigh_nodes(set<int> & connected_nodes,NodeFracturesMap node_fract_map,int node_index){
	for (HydrauPropMap::iterator it1=node_fract_map[node_index].begin();it1!=node_fract_map[node_index].end();it1++){
		if (connected_nodes.find(it1->first)==connected_nodes.end()){
			connected_nodes.insert(it1->first);
			add_neigh_nodes(connected_nodes,node_fract_map,it1->first);
		}
	}
}


/*
void VerticalUpdateFracture(NetworkMeshes& net_mesh,vector<FractureMesh>& new_meshes,int i,double y1,double y2){
	// variables
	FluxPoint2D p1,p2;
	int index1,index2;
	BordersMap::iterator it1;NodesMap::iterator it2;
	// update fracture description
	index1=net_mesh.meshes[i].p_ori.index;index2=net_mesh.meshes[i].p_tar.index;
	p1=FluxPoint2D(net_mesh.meshes[i].p_ori.p.x(),y1,index1);
	p2=FluxPoint2D(net_mesh.meshes[i].p_tar.p.x(),y2,index2);
	new_meshes.push_back(FractureMesh(p1,p2,net_mesh.meshes[i].velocity,net_mesh.meshes[i].aperture,net_mesh.meshes[i].mesh_index,net_mesh.meshes[i].neigh_meshes));
	// update border map
	it1=net_mesh.border_map.find(index1);
	if (it1!=net_mesh.border_map.end()){net_mesh.border_map[index1]=it1->second;}
	it1=net_mesh.border_map.find(index2);
	if (it1!=net_mesh.border_map.end()){net_mesh.border_map[index2]=it1->second;}
	// update node map
	it2=net_mesh.nodes_map.find(index1);
	if (it2!=net_mesh.nodes_map.end()){net_mesh.nodes_map[index1]=p1.p;}
	it2=net_mesh.nodes_map.find(index2);
	if (it2!=net_mesh.nodes_map.end()){net_mesh.nodes_map[index2]=p2.p;}
}
*/

NetworkMeshes Vertical_Translation(NetworkMeshes net_mesh,double H){
	// 1. Variables
	double x1,x2,y1,y2,min_y=net_mesh.domain.min_pt.j,max_y=net_mesh.domain.max_pt.j;
	if (H<0){cout << "WARNING1 in Vertical_Translation (NetworkMeshes.cpp): case not implemented" << endl;}
	FractureMesh frac_mesh;
	// 2. New network initialization
	NetworkMeshes new_net_mesh;new_net_mesh.domain=net_mesh.domain;
	new_net_mesh.max_fract_spacing=net_mesh.max_fract_spacing;
	// 3. Loop on the initial fractures for translation and addition to the new network
	for (int i=0;i<net_mesh.meshes.size();i++){
		// 3.1. Translation
		// first extremity
		x1=net_mesh.meshes[i].p_ori.p.x();
		y1=net_mesh.meshes[i].p_ori.p.y()+H;
		// second extremity
		x2=net_mesh.meshes[i].p_tar.p.x();
		y2=net_mesh.meshes[i].p_tar.p.y()+H;
		// 3.2. Addition to the new network
		// translated fracture is still in the domain
		if (y1<max_y&&y2<max_y){
			frac_mesh=FractureMesh(FluxPoint2D(x1,y1),FluxPoint2D(x2,y2),net_mesh.meshes[i].aperture,net_mesh.meshes[i].fracture_index);
			AddFracture(frac_mesh,new_net_mesh);
		}
		// fracture totally out of the domain
		else if (y1>max_y&&y2>max_y){
			y1=y1-max_y+min_y;
			y2=y2-max_y+min_y;
			frac_mesh=FractureMesh(FluxPoint2D(x1,y1),FluxPoint2D(x2,y2),net_mesh.meshes[i].aperture,net_mesh.meshes[i].fracture_index);
			AddFracture(frac_mesh,new_net_mesh);
		}
		// piece of fracture out of the domain, case 1
		else if (y1<max_y&&y2>max_y){
			if (fabs(x1-x2)>EPSILON){cout << "WARNING1 in Vertical_Translation (NetworkMeshes.cpp): case not implemented" << endl;}
			// first segment
			frac_mesh=FractureMesh(FluxPoint2D(x1,y1),FluxPoint2D(x1,max_y),net_mesh.meshes[i].aperture,net_mesh.meshes[i].fracture_index);
			AddFracture(frac_mesh,new_net_mesh);
			// second segment
			y2=y2-max_y+min_y;
			frac_mesh=FractureMesh(FluxPoint2D(x1,min_y),FluxPoint2D(x1,y2),net_mesh.meshes[i].aperture,net_mesh.meshes[i].fracture_index);
			AddFracture(frac_mesh,new_net_mesh);
		}
		// piece of fracture out of the domain, case 2
		else if (y1>max_y&&y2<max_y){
			if (fabs(x1-x2)>EPSILON){cout << "WARNING1 in Vertical_Translation (NetworkMeshes.cpp): case not implemented" << endl;}
			// first segment
			frac_mesh=FractureMesh(FluxPoint2D(x1,y2),FluxPoint2D(x1,max_y),net_mesh.meshes[i].aperture,net_mesh.meshes[i].fracture_index);
			AddFracture(frac_mesh,new_net_mesh);
			// second segment
			y1=y1-max_y+min_y;
			frac_mesh=FractureMesh(FluxPoint2D(x1,min_y),FluxPoint2D(x1,y1),net_mesh.meshes[i].aperture,net_mesh.meshes[i].fracture_index);
			AddFracture(frac_mesh,new_net_mesh);
		}
		else{cout << "WARNING2 in Vertical_Translation (NetworkMeshes.cpp): case not implemented" << endl;}
	}
	return new_net_mesh;
}

void FractureMeshNumbering(NetworkMeshes& net_mesh){
	// 1. Define the number of each fracture mesh
	int cpt_mesh=-1;
	for (int i=0;i<net_mesh.meshes.size();i++){
		cpt_mesh++;
		net_mesh.meshes[i].mesh_index=cpt_mesh;
		if (net_mesh.meshes[i].original_mesh_index<0){
			net_mesh.meshes[i].original_mesh_index=cpt_mesh;
		}
	}
	net_mesh.cpt_mesh=cpt_mesh;
	// 2. Define the list of neighbors of each mesh
	int index1,index2;
	for (int i=0;i<net_mesh.meshes.size();i++){
		index1=net_mesh.meshes[i].p_ori.index;
		index2=net_mesh.meshes[i].p_tar.index;
		for (int j=0;j<net_mesh.meshes.size();j++){
			if (i!=j){
				if (net_mesh.meshes[j].p_ori.index==index1||net_mesh.meshes[j].p_tar.index==index1||net_mesh.meshes[j].p_ori.index==index2||net_mesh.meshes[j].p_tar.index==index2){
					net_mesh.meshes[i].neigh_meshes.push_back(net_mesh.meshes[j].mesh_index);
				}
			}
		}
	}
}

void UpdateCptFract(NetworkMeshes& net_mesh){
	set<int> fracture_indices;
	// 1. Redefine fracture_indices depending on the meshes that disappeared when computing the backbone
	for (int i=0;i<net_mesh.meshes.size();i++){
		fracture_indices.insert(net_mesh.meshes[i].fracture_index);
	}
	// 2. Redefine fracture_indices depending on the fractures that disappeared when computing the backbone
	set<int>::iterator it_set=fracture_indices.end();it_set--;
	if (fracture_indices.size()-1!=*it_set){
		vector<FractureMesh> new_meshes=net_mesh.meshes;
		int cpt_fract=0;
		for (set<int>::iterator it_fract=fracture_indices.begin();it_fract!=fracture_indices.end();it_fract++){
			for (int i=0;i<net_mesh.meshes.size();i++){
				if (net_mesh.meshes[i].fracture_index==*it_fract){new_meshes[i].fracture_index=cpt_fract;}
			}
			cpt_fract++;
		}
		net_mesh.meshes=new_meshes;
	}
	net_mesh.cpt_fract=fracture_indices.size()-1;
	net_mesh.fract_indices=fracture_indices;
}

/********************************* FUNCTIONS TO COMPUTE THE FLOW VELOCITY ***********************/
double VelocityComputation(double aperture,double h1,double h2,double length){
	double transmissivity=ReturnTransmissivity(aperture);
	return transmissivity/aperture*(h1-h2)/length;
}

// Compute flow and return the fracture network composed of fracture with flow
NetworkMeshes FlowComputation(NetworkMeshes net_mesh,BoundaryConditionsDef bc_map_def){
	NetworkMeshes flow_net_mesh=net_mesh;
	// 1. Boundary condition definition
	BoundaryConditionsDFN bc_map=ReturnBoundCondDFN(net_mesh.domain,bc_map_def,net_mesh.border_map,net_mesh.nodes_map);
	// 2. Pressure computation
	ublas_vector hydraulic_head=DFNComputationClassic(net_mesh,bc_map);

	//print_vector(hydraulic_head,hydraulic_head.size());

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
