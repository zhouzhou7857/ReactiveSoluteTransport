/*
 * Projection.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */


#include "Transport.h"
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

// Projection of a point within the domain
void Transport::M_project_in_system(pointcpp<double> & M_out, const double & EPS){
	// if M is out the system, project it inside the system
	pointcpp<double> pmin,pmax;
	domain.get_extremities(pmin,pmax);
	if(M_out.i < pmin.i-EPS) M_out.i = pmin.i + EPS;
	if(M_out.i > pmax.i+EPS) M_out.i = pmax.i - EPS;
	if(M_out.j < pmin.j-EPS) M_out.j = pmin.j + EPS;
	if(M_out.j > pmax.j+EPS) M_out.j = pmax.j - EPS;
}

// Projection of a point within the domain
void Transport::M_project_in_system(CgalPoint2D & M_out, const double & EPS){
	// conversion from CgalPoint2D to pointcpp
	pointcpp<double> point(M_out);
	// projection in the system
	M_project_in_system(point,EPS);
	// conversion from pointcpp to CgalPoint2D
	M_out = point.CgalPoint();
}

/**
*returns the projection through the domain boundaries (consider periodic boundaries)
*when there is no projection on seg_ortho1 -> we use the projection by periodicity on seg_ortho2
*/
bool Transport::Projection_Periodic_Boundaries(CgalPoint2D init_point,Segment2D seg_ortho1,Segment2D seg_ortho2,CgalPoint2D M_end1,CgalPoint2D M_end2,CgalPoint2D & proj_point,int & proj_mesh,double & proj_dist){
	// 0. Parameters
	bool success = false;
	vector<FractureMesh> meshes = this->net_mesh.meshes;
	// 1. Verification of the closest projection
	if (Closest_Orthogonal_Projection(init_point,seg_ortho1,proj_point,proj_mesh,proj_dist)){
		std::cout<<"\n PB IN Projection_Periodic_Boundaries (Projection.cpp): closest projection exists";
	}
	// 2. Determination of the projection through the periodic boundary
	else{
		// 2.1. Determination of the projection point -> the further from the initial point (closest to the second boundary)
		CgalPoint2D current_point;
		double max_dist = -1;
		int meshes_size = meshes.size();
		for (int i = 0; i<meshes_size; i++){
			if ((!meshes[i].Is_Collinear(init_point))&&(meshes[i].Is_Intersected(seg_ortho2, current_point))&&(!identic(init_point,current_point,EPSILON))&&(distance_2D(init_point,current_point)>max_dist)){
				proj_mesh = meshes[i].mesh_index;
				proj_point = current_point;
				max_dist = distance_2D(init_point,proj_point);
				success = true;
			}
		}
		// 2.2. Determination of the distance between the initial point to the projected point through the domain periodic boundaries
		// 2.2.1. Distance from the initial point to the closest boundary
		M_project_in_system(M_end1,EPSILON);
		double dist1 = distance_2D(init_point,M_end1);
		// 2.2.2. Distance from the second boundary to the projected point
		M_project_in_system(M_end2,EPSILON);
		double dist2 = distance_2D(M_end2,proj_point);
		// 2.2.3. total distance
		proj_dist = dist1+dist2;
	}
	return success;
}

/**returns the definition (point, mesh index and distance) of the closest projection of M orthogonal to seg_ortho*/
bool Transport::Closest_Orthogonal_Projection(CgalPoint2D init_point,Segment2D seg_ortho,CgalPoint2D & proj_point,int & proj_mesh,double & proj_dist){
	// 0. Parameters
	bool succes = false;
	vector<FractureMesh> meshes = this->net_mesh.meshes;
	int meshes_size = meshes.size();
	// 1. loop on the meshes to determine the closest intersection
	proj_dist = 1e9;
	CgalPoint2D current_point;
	for (int i = 0; i<meshes_size; i++){
		if ((!meshes[i].Is_Collinear(init_point))&&(meshes[i].Is_Intersected(seg_ortho, current_point))&&(!identic(init_point,current_point,EPSILON))&&(distance_2D(init_point,current_point)<proj_dist)){
			proj_mesh = meshes[i].mesh_index;
			proj_point = current_point;
			proj_dist = distance_2D(init_point,proj_point);
			succes = true;
		}
	}
	return succes;
}

double Transport::ReturnPureDiffusionTime(double distance){
	double u_simu = rng_tracker.uniform(0,1);
	double value=distance/(2*sqrt(phys_param.Dm)*boost::math::erf_inv(1-u_simu));
	return pow(value,2.);
}

/*
 * For a particle reaching going from Minit (M1) to Mproj (Mborder) with Mproj on a domain border, returns its new position, new fracture mesh and required time to reach it
 * * considering a reflecting condition on the top and bottom border and displacement by diffusion in the matrix
 * * returns true if the particle reaches the right border of the domain, and false otherwise
 */
bool Transport::ReflectionOnBorder(pointcpp<double>& Minit,pointcpp<double>& Mproj,int & transfer_mesh,double & transfer_time){
	// 0. Check that the extremity of seg_proj is on one of the domain border
	CgalPoint2D M_border=Mproj.CgalPoint();
	string reflec_border=domain.ReturnBorder(M_border);
	if (reflec_border==NO_BORDER||reflec_border==RIGHT_BORDER){cout << "PB in ReflectionOnBorder (Projection.cpp): a border projection is not defined" << endl;exit(1);}
	// 1. Compute the symmetric segment and line
	CgalPoint2D M1=Minit.CgalPoint(),M1_sym=domain.ReturnSymmetric(M1,M_border);
	Segment2D seg_sym(M_border,M1_sym);CgalLine2D line_sym=seg_sym.supporting_line();
	// 2. Determine the intersections of this symmetric line
	CgalPoint2D M_inter1,M_inter2,M_inter_new;
	if (!domain.IntersectionBorders(line_sym,M_inter1,M_inter2)||!define(M_inter1)||!define(M_inter2)){
		// new search of the two intersections (due to numerical issue with CGAL)
		if (!domain.IntersectionBordersExtended(line_sym,M_inter1,M_inter2)||!define(M_inter1)||!define(M_inter2)){
			std::cout<<"\n PB IN ReflectionOnBorder (Projection.cpp): intersections with domain borders not found" << endl;exit(1);
		}
	}
	if (identic(M_inter1,M_border, EPSILON)){M_inter_new=M_inter2;}
	else if (identic(M_inter2,M_border, EPSILON)){M_inter_new=M_inter1;}
	else{std::cout<<"\n PB IN ReflectionOnBorder (Projection.cpp): intersections not well defined" << endl;exit(1);}
	// 3. Look for the next projection (the new Mproj)
	Minit=Mproj;
	// 3.1. If a fracture intersects the new direction
	Segment2D seg_inter(M_border,M_inter_new);CgalPoint2D inter_fract;double transfer_dist;
	if (Closest_Orthogonal_Projection(M_border,seg_inter,inter_fract,transfer_mesh,transfer_dist)){
		Mproj=pointcpp<double>(inter_fract);
		transfer_time=ReturnPureDiffusionTime(transfer_dist);
		return false;
	}
	// 3.2. Else, the new position is on one of the domain border
	else{
		Mproj=M_inter_new;
		transfer_dist=distance_2D(M_border,M_inter_new);
		transfer_time=ReturnPureDiffusionTime(transfer_dist);
		transfer_mesh=-1;
		// If the new position is on the right border
		if (domain.ReturnBorder(M_inter_new)==RIGHT_BORDER){return true;}
		//  If the new position is on another border
		else if (domain.ReturnBorder(M_inter_new)!=NO_BORDER){return false;}
	}
	std::cout<<"\n PB IN ReflectionOnBorder (Projection.cpp): case not implemented" << endl;exit(1);
	return false;
}

/**
* returns the two closest projections (on the fractures located below and above)  of M1 orthogonal to the segment (M1,M2)
* map<projected point,pair<mesh index,distance from initial point to projected point>>
*/
//Projection_Map Transport::Orthogonal_Projection_M(CgalPoint2D M1, CgalPoint2D M2){
bool Transport::Orthogonal_Projection_M(CgalPoint2D M1, CgalPoint2D M2,Projection_Map& proj_points){
	//0. PARAMETERS
	//Projection_Map proj_points;
	CgalPoint2D inter;
	//1. CONSTRUCTION OF THE SEGMENTS (ONE BY SIDE) ORTHOGONAL TO SEG AND STARTING FROM THE MIDDLE OF THE SEG
	Segment2D seg(M1,M2);CgalLine2D line=seg.supporting_line();//1.1. Definition of the initial segment and its supporting line
	CgalLine2D line_ortho=line.perpendicular(M1);CgalPoint2D M_end1,M_end2;//1.2. Definition of the perpendicular line and segments
	if (!domain.IntersectionBorders(line_ortho,M_end1,M_end2)||!define(M_end1)||!define(M_end2)){
		// new search of the two intersections (due to numerical issue with CGAL)
		if (!domain.IntersectionBordersExtended(line_ortho,M_end1,M_end2)||!define(M_end1)||!define(M_end2)){
			//std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): intersections with domain borders not found";// << endl;
			/*if (!domain.IntersectionBordersExtended(line_ortho,M_end1,M_end2)){cout << "Problem 1" << endl;}
			if (!define(M_end1)){cout << "Problem 2" << endl;}
			if (!define(M_end2)){cout << "Problem 3" << endl;}
			cout << "M1 = ";print(M1);
                        cout << "M2 = ";print(M2);*/
			//exit(1);
			return false;
		}
	}
	Segment2D seg_ortho1(M1,M_end1),seg_ortho2(M1,M_end2);
	//2. DETERMINATION OF THE INTERSECTIONS WITH THE FRACTURE NETWORK
	double min_dist = 1e9;int mesh_index = -1;
	//2.1. determination of the projection : the closest fracture or boundary (considering reflecting boundaries)
	//2.1.1. first direction: from M to Mend1
	/*if (!Closest_Orthogonal_Projection(M1,seg_ortho1,inter,mesh_index,min_dist)){
		if (domain.ReturnBorder(seg_ortho1.source())!=NO_BORDER&&domain.ReturnBorder(seg_ortho1.target())==NO_BORDER){inter=seg_ortho1.source();}
		else if (domain.ReturnBorder(seg_ortho1.source())==NO_BORDER&&domain.ReturnBorder(seg_ortho1.target())!=NO_BORDER){inter=seg_ortho1.target();}
		else{std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;exit(1);}
		mesh_index = -1;min_dist = distance_2D(M1,inter);// update values
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);*/
	if (Closest_Orthogonal_Projection(M1,seg_ortho1,inter,mesh_index,min_dist)){
		//proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);
		proj_points[pointcpp<double>(inter)] = std::make_pair(mesh_index,min_dist);
	}
	//2.1.2. second direction: from M to Mend2
	/*if (!Closest_Orthogonal_Projection(M1,seg_ortho2,inter,mesh_index,min_dist)){
		if (domain.ReturnBorder(seg_ortho2.source())!=NO_BORDER&&domain.ReturnBorder(seg_ortho2.target())==NO_BORDER){inter=seg_ortho2.source();}
		else if (domain.ReturnBorder(seg_ortho2.source())==NO_BORDER&&domain.ReturnBorder(seg_ortho2.target())!=NO_BORDER){inter=seg_ortho2.target();}
		else{std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;exit(1);}
		mesh_index = -1;min_dist = distance_2D(M1,inter);// update values
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);*/
	if (Closest_Orthogonal_Projection(M1,seg_ortho2,inter,mesh_index,min_dist)){
		//proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);
                proj_points[pointcpp<double>(inter)] = std::make_pair(mesh_index,min_dist);
	}
	// if there is not 2 intersections -> problem
	//if (proj_points.size()!=2){cout << "\n PB IN Orthogonal_Projection_M (Projection.cpp) : there is not two intersections" << endl;exit(1);}
	//return proj_points;
	return true;
}


/*Projection_Map Transport::Orthogonal_Projection_M(pointcpp<double> M, pointcpp<double> M2){
	//0. PARAMETERS
	Projection_Map proj_points;CgalPoint2D inter;
	double L = max(this->domain.domain_size_x(),this->domain.domain_size_y());
	//1. CONSTRUCTION OF THE SEGMENTS (ONE BY SIDE) ORTHOGONAL TO SEG AND STARTING FROM THE MIDDLE OF THE SEG
	//1.1. middle and angle of the initial segment
	CgalPoint2D M1 = M.CgalPoint(), M2_ = M2.CgalPoint();
	Segment2D seg(M1,M2_);
	double alpha = seg.get_orientation();	//angle of the current mesh
	//1.2. determination of the extremities of the orthogonal segments
	alpha += PI/2;
	CgalPoint2D M_end1 = CgalPoint2D(M1.x()+2*L*std::cos(alpha),M1.y()+2*L*std::sin(alpha));
	CgalPoint2D M_end2 = CgalPoint2D(M1.x()-2*L*std::cos(alpha),M1.y()-2*L*std::sin(alpha));
	Segment2D seg_ortho1(M1,M_end1),seg_ortho2(M1,M_end2);
	//2. DETERMINATION OF THE INTERSECTIONS WITH THE FRACTURE NETWORK
	//2.0. Parameters to describe the selected projection
	double min_dist = 1e9;int mesh_index = -1;
	//2.1. determination of the projection : the closest one or the one through periodic boundary
	//2.1.1. first direction: from M to Mend1
	// Considering an infinite matrix around the studied domain
	/*if (Closest_Orthogonal_Projection(M1,seg_ortho1,inter,mesh_index,min_dist)){
		proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);
	}*/
	// Considering periodic boundary conditions on the top and bottom border
	/*if (!Closest_Orthogonal_Projection(M1,seg_ortho1,inter,mesh_index,min_dist)){
		if (!Projection_Periodic_Boundaries(M1,seg_ortho1,seg_ortho2,M_end1,M_end2,inter,mesh_index,min_dist)){
			std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;
			exit(1);
		}
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);*/
	// Considering reflecting boundaries at the domain borders
/*	if (!Closest_Orthogonal_Projection(M1,seg_ortho1,inter,mesh_index,min_dist)){
		CgalPoint2D inter_border2;
		if (!domain.IntersectionBorders(seg_ortho1,inter,inter_border2)){
			std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;
			exit(1);
		}
		// check values for the intersection with the domain borders
		if (inter_border2.x()!=NOT_DEFINED||inter_border2.y()!=NOT_DEFINED){
			std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): a second border intersection exists" << endl;
		}
		// update values
		mesh_index = -1;min_dist = distance_2D(M1,inter);
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);

	//2.1.2. second direction: from M to Mend2
	if (Closest_Orthogonal_Projection(M1,seg_ortho2,inter,mesh_index,min_dist)){
		proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);
	}
	/*if (!Closest_Orthogonal_Projection(M1,seg_ortho2,inter,mesh_index,min_dist)){
		if (!Projection_Periodic_Boundaries(M1,seg_ortho2,seg_ortho1,M_end2,M_end1,inter,mesh_index,min_dist)){
			std::cout<<"\n PB IN Orthogonal_Projection_M (Projection.cpp): orthogonal projection not determined" << endl;
		}
	}
	proj_points[pointcpp<double>(inter)] = std::make_pair<int,double>(mesh_index,min_dist);*/
	// if there is not 2 intersections -> problem
	/*if (proj_points.size()!=2){
		/*cout << "\n PB IN Orthogonal_Projection_M (Projection.cpp) : method implemented only for 2 intersections" << endl;
		cout << "Number of intersections is " << proj_points.size() << endl;*/
	//}
/*	return proj_points;
}*/
