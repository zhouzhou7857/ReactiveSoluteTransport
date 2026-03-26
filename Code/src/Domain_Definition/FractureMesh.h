/*
 * FractureMesh.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef FRACTUREMESH_H_
#define FRACTUREMESH_H_

#include "../Utilitaries/FluxPoint2D.h"
#include "../Utilitaries/Pointcpp.h"
#include "../Utilitaries/RandomNumber.h"
#include "../Utilitaries/Segment.h"
//#include "../Utilitaries/Constantes.h"
#include "../Utilitaries/Structures.h"


// definition of a mesh
class FractureMesh{
public:
	FluxPoint2D p_ori; /**< origin of the mesh. */
	FluxPoint2D p_tar;	/**< target of the mesh. */
	double velocity;	/**< the velocty in the mesh */
	double aperture;	// mesh aperture
	int mesh_index; /**< the index of the mesh .*/
	int original_mesh_index; /**< stable mesh index from initial numbering */
	std::vector<int> neigh_meshes;	/**< the neighboring meshes.*/
	double previous_beta;
	int fracture_index;
public:
	FractureMesh(){};
	virtual ~FractureMesh(){};
	FractureMesh(FluxPoint2D,FluxPoint2D,double,double,int,std::vector<int>,int);
	FractureMesh(FluxPoint2D,FluxPoint2D,double,double,int);
	FractureMesh(FluxPoint2D,FluxPoint2D,double,int);
	FractureMesh(FluxPoint2D,FluxPoint2D,double,int,int);
	FractureMesh(Segment2D,double,int);
	double distance_from_M_to_extremity(pointcpp<double>, pointcpp<double> &);
	double Mesh_Time_Advection(pointcpp<double>,pointcpp<double>);
	bool Get_Total_Time_From_Advec_Time1D(double&,double,double,double,RngStream_a &,double);
	//bool Get_Total_Time_From_Advec_Time2D(double&,double,double,double,RngStream_a &);
	// VRG Dec 2016 implementing Method H
	//bool Get_Total_Time_From_Advec_Time2D_Method_H(double&,double,double,double,RngStream_a &,double,double);
	// VRG March 2017 implementing iteration algorithm Method H iAiiA with patches to address lost particles
	//bool Get_Total_Time_From_Advec_Time2D_HT2(double&,double,double,double,RngStream_a &,double,double);
	bool Get_Total_Time_From_Advec_Time2D_Conv(double&,double,double,double,RngStream_a &,double,double,double,double,Declare_CDF_Map&,int);
	//double ReturnDiffusionTime_final(double,double,double,double, double, double, double, double, double, double);
	double ReturnDiffusionTime_HT2(double,double,double,double, double, double, double, double, double, double);
	double ReturnDiffusionTime_bisection(double,double,double,double, double, double, double, double, double, double);
	//double ReturnDiffusionTime_bisection_b(double,double,double,double, double, double, double, double, double, double);
	double ReturnDiffusionTime_Con (double&,double,double,double,RngStream_a &,double,double,double,double,Declare_CDF_Map&,int);	// VRG Nov 2017
	double ReturnDiffusionTime_delta(std::vector<double>, std::vector<double>, double);
	//double ReturnDiffusionTime_delta(double *Time_Dist, std::vector<double>, double);
	std::vector<double> ReturnPMF_Delta(std::vector<double>, std::vector<double> );	// VRG Nov 2017
	std::vector<double> ReturnCDF_Delta(std::vector<double>); // VRG March 2018
	std::vector<double> Shave_PMF(std::vector<double> ,int ); // VRG March 2018
	//std::vector<double>  ReturnCDF_Delta(double,double,double,double,double,double,double,double,double,double,double,double *Time_Dist); // VRG Nov 2017
	std::vector<double>  ReturnCDF_Delta(double,double,double,double,double,double,double,double,double,double,double,std::vector<double>); // VRG Nov 2017
	bool Get_Total_Time_From_HT1(double,double&,double,double,double,double);	// VRG Nov 2017



	//bool Get_Total_Time_From_Advec_Time(double&,double,double,double,RngStream_a &,double,double,double,Declare_CDF_Map&,int);

	double Advection_Time_Computation(pointcpp<double>,const Projection_Map&,const double,const double,const double,const double);
	int Is_Collinear(CgalPoint2D);
	int Is_Intersected(Segment2D&,CgalPoint2D&);
	pointcpp<double> velocity_interpolation(const pointcpp<double>&);
	pointcpp<double> Mesh_Scale_Advection(pointcpp<double>,double,bool opposite=false);
	Segment2D ReturnSegment();
	double ReturnLength();
	double ReturnConductance();
	// modif delphine - 2016/12/09
	//double Temp_fnc1(double,double,double,double,double,double,double,double,double,double);
	//double Temp_fnc2(double,double,double,double,double,double,double,double,double);
	double Temp_fncb(double,double,double,double,double,double,double,double,double,double);
	double Temp_fncU(double,double,double,double,double,double,double,double,double,double);
	double Temp_fnc_max(double,double, double,double,double,double,double,double);
	double Temp_fnc_c(double,double, double,double,double,double,double,double,double);
};

double sign(double);


#endif /* FRACTUREMESH_H_ */
