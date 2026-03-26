/*
 * FractureMesh.cpp
 *		HEAT
 *  	Created on: August 31, 2015   & Updated for 2D diffusion October 2016
 *      Author: viktoria and delphine
 */


#include "FractureMesh.h"
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <limits>
#include <fstream>


using namespace std;


FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double velocity_,double aperture_,int mesh_index_,std::vector<int> neigh_meshes_,int fracture_index_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	velocity = velocity_;
	aperture = aperture_;
	mesh_index = mesh_index_;
	original_mesh_index = mesh_index_;
	neigh_meshes = neigh_meshes_;
	fracture_index=fracture_index_;
}

FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double velocity_,double aperture_,int mesh_index_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	velocity = velocity_;
	aperture = aperture_;
	mesh_index = mesh_index_;
	original_mesh_index = mesh_index_;
}

FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double aperture_,int mesh_index_, int fracture_index_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	aperture = aperture_;
	mesh_index = mesh_index_;
	original_mesh_index = mesh_index_;
	fracture_index=fracture_index_;
}

FractureMesh::FractureMesh(FluxPoint2D p_ori_,FluxPoint2D p_tar_,double aperture_,int fracture_index_){
	p_ori = p_ori_;
	p_tar = p_tar_;
	aperture = aperture_;
	original_mesh_index = -1;
	fracture_index=fracture_index_;
}

FractureMesh::FractureMesh(Segment2D seg,double aperture_,int fracture_index_){
	p_ori=FluxPoint2D(seg.source());
	p_tar=FluxPoint2D(seg.target());
	aperture = aperture_;
	original_mesh_index = -1;
	fracture_index=fracture_index_;
}


// return the distance from the current position M to the extremity for a displacement with the velocity
double FractureMesh::distance_from_M_to_extremity(pointcpp<double> M, pointcpp<double> & M_out){
	double distance = 0;
	CgalPoint2D M_(M.i,M.j);
	if (velocity>0){
		distance = sqrt(CGAL::squared_distance(M_,p_tar.p));
		M_out.define((double)p_tar.p.x(),(double)p_tar.p.y());
	}
	else if (velocity<0){
		distance = sqrt(CGAL::squared_distance(M_,p_ori.p));
		M_out.define((double)(p_ori.p.x()),(double)(p_ori.p.y()));
	}
	return distance;
}

// return the time required to travel by advection from M_in to M_out in the current mesh
double FractureMesh::Mesh_Time_Advection(pointcpp<double> M_in, pointcpp<double> M_out){

	if (velocity == 0){
		cout << "WARNING in Mesh_Time_Advection (FractureMesh.cpp): velocity null" << endl;
	}

	CgalPoint2D M_in_(M_in.i,M_in.j), M_out_(M_out.i,M_out.j);
	double distance = sqrt(CGAL::squared_distance(M_in_,M_out_));
	return distance/fabs(velocity);
}

double sign(double x){
	if (x<0){return -1;}
	else if (x==0){return 0;}
	else if (x>0){return 1;}
	cout << "WARNING in sign (FractureMesh.cpp): x is not defined" << endl;
	return -1;
}

// Returns the maximum temperature that can be reached at a point L with a given maximum arrival time
double FractureMesh::Temp_fnc_max(double r_0,double r,double r_1,double Dm,double velocity,double advection_time,double arr_time_max,double L){
	double advection_time_L=L/fabs(velocity);  // Currently the velocity is the same for the whole fracture not just segment so we can use this approach
	double value1=-(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)/(4*Dm*arr_time_max);
	double value2=boost::math::expint(value1);
	double result=1./(2*r_0)*(sign(fabs(velocity)*advection_time_L)-(r/PI)*value2-erf((r*fabs(velocity)*advection_time_L)/(2*sqrt(r_0*Dm*arr_time_max))))+r_1/(2*PI*pow(r_0,1.5)*advection_time_L)*
			exp(-(r*r*Dm*arr_time_max)/(4*r_0*(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)))*(2*advection_time_L/r_1-sqrt(PI*arr_time_max*Dm)/(r*r_0*fabs(velocity)*advection_time_L)-1/(r_0*r_0));
	if(result<0){ cout << "WARNING in Temp_fnc_max (FractureMesh.cpp):  is negative" << endl; }

	//if(fabs(result-0.846731)<1e2){ cout << "T_max" << result << "advection_time" <<advection_time << "velocity" << velocity << "L" << L << "arr_time_max" << arr_time_max << endl; }
	return result;

}

// March 30 2017 VRG: added Temp_fncU so we can find the root of f(t)/f(t_max)-U=0 for the purpose of implementing Method H iAiiA.
double FractureMesh::Temp_fncU(double r_0,double r,double r_1,double velocity,double Dm,double beta,double time,double L, double U, double Temp_max){
	double advection_time_L=L/fabs(velocity);  // // We cant use advection_time becaus that is for the segment we need the advection time all the way from the orgin hence advection_time_L
	double value1=-(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)/(4*Dm*(time+advection_time_L));
	double value2=boost::math::expint(value1);
	double result=1./(2*r_0*Temp_max)*(sign(fabs(velocity)*advection_time_L)-(r/PI)*value2-erf((r*fabs(velocity)*advection_time_L)/(2*sqrt(r_0*Dm*(time+advection_time_L)))))+r_1/(2*PI*pow(r_0,1.5)*advection_time_L)*
			exp(-(r*r*Dm*(time+advection_time_L))/(4*r_0*(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)))*(2*advection_time_L/r_1-sqrt(PI*(time+advection_time_L)*Dm)/(r*r_0*fabs(velocity)*advection_time_L)-1/(r_0*r_0))-U;

	return result;
}

// March 30 2017 VRG: added Temp_fncb so we can find the root of f(t)-beta=0 for the purpose of implementing Method H iAiiA.
double FractureMesh::Temp_fncb(double r_0,double r,double r_1,double velocity,double Dm,double beta,double time,double L, double U, double Temp_max){
	double advection_time_L=L/fabs(velocity);  // // We can't use advection_time because that is for the segment we need the advection time all the way from the origin hence advection_time_L
	double value1=-(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)/(4*Dm*(time+advection_time_L));
	double value2=boost::math::expint(value1);
    double result=1./(2*r_0)*(sign(fabs(velocity)*advection_time_L)-(r/PI)*value2-erf((r*fabs(velocity)*advection_time_L)/(2*sqrt(r_0*Dm*(time+advection_time_L)))))+r_1/(2*PI*pow(r_0,1.5)*advection_time_L)*
			exp(-(r*r*Dm*(time+advection_time_L))/(4*r_0*(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)))*(2*advection_time_L/r_1-sqrt(PI*(time+advection_time_L)*Dm)/(r*r_0*fabs(velocity)*advection_time_L)-1/(r_0*r_0))-beta;
	return result;
}

// Nov 2017 the function used for HT2 convolution approach.
double FractureMesh::Temp_fnc_c(double r_0,double r,double r_1,double Dm,double velocity,double advection_time,double arr_time_max,double L,double Temp_max){
	double advection_time_L=L/fabs(velocity);  // Currently the velocity is the same for the whole fracture not just segment so we can use this approach
	double value1=-(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)/(4*Dm*arr_time_max);
	double value2=boost::math::expint(value1);
	double result=1./(2*r_0*Temp_max)*(sign(fabs(velocity)*advection_time_L)-(r/PI)*value2-erf((r*fabs(velocity)*advection_time_L)/(2*sqrt(r_0*Dm*arr_time_max))))+r_1/(2*PI*pow(r_0,1.5)*advection_time_L)*
			exp(-(r*r*Dm*arr_time_max)/(4*r_0*(fabs(velocity)*advection_time_L)*(fabs(velocity)*advection_time_L)))*(2*advection_time_L/r_1-sqrt(PI*arr_time_max*Dm)/(r*r_0*fabs(velocity)*advection_time_L)-1/(r_0*r_0));

	return result;

}


// VRG Dec 2016 - This is also part of the implementation of Method H.
// VRG modified March 2017 -- the bisection algorithm now uses definition f(t)/f(t_max)-U=0 in its root finding trough the function Temp_fncU.
double FractureMesh::ReturnDiffusionTime_bisection(double r_0,double r,double r_1,double velocity, double beta, double L, double Dm, double arr_time_max, double U, double Temp_max){
	    double t_max=MAX_TIME;// ;  // This is our observation time of the system e.g. 50 years
	   double eps=std::numeric_limits<double>::epsilon(); double tau=sqrt(eps); // squaroot of the machine epsilon is the tolerance for error
	    double tmin=0, tmid, tmax, count=0;

	    tmax=U*t_max*1.1;
	    double time =tmin;
	    double flag=0;
	    while (flag<1){
	        tmid=(tmin+tmax)/2;         //here we have a bisection method to assist if the value of x_k+1 is out of bound
	        // here the second control addresses the case where time is practically zero
	        if (abs(Temp_fncU(r_0,r,r_1,fabs(velocity),Dm,beta,tmid,L,U,Temp_max))< eps) {
        								flag=2;
        				        		if ( t_max < tmid){
        				        	        cout << " Warning (FractureMesh.cpp): Control 1 in bisection case a fails ........A" << endl;
        				        		}
	        }
	        if (abs((tmax-tmin)/tmax) < eps && tmid<t_max ){
	        	       flag=4;
	        		if ( t_max < tmid){
	        	        cout << " Warning (FractureMesh.cpp): tmid= " << tmid << endl;
	        	        cout << " Warning (FractureMesh.cpp): U= " << U<< endl;
	        	        cout << " Warning (FractureMesh.cpp): Control 1 in bisection case b fails ........B" << endl;
	        		}
	        }
	       count=count+1;
	       if (count > 5000){
	    	   	   	flag=3;
    		        cout << " Warning (FractureMesh.cpp):The root finding algorithm failed the bisection method did not address the numerical error see ReturnDiffusionTime_bisection" << endl;
	       }
	       if (sign(Temp_fncU(r_0,r,r_1,fabs(velocity),Dm,beta,tmid,L,U,Temp_max)) == sign(Temp_fncU(r_0,r,r_1,fabs(velocity),Dm,beta,tmin,L,U,Temp_max))){
	    	   	   tmin=tmid;}
	    	   else{tmax=tmid;}
	    }
	    time=tmid;
	    return time;
}

// VRG March 2017 - implementing Method H vs iAiiA root finding algorithm with the patch that addresses lost particles.
double FractureMesh::ReturnDiffusionTime_HT2(double r_0,double r,double r_1,double velocity, double beta, double L, double Dm, double arr_time_max, double U, double Temp_max){

	double t_max=MAX_TIME;// ;  // This is our observation time of the system e.g. 50 years

	    double xmin=1;
		double xmax=8.6e+127;
		double time=xmin;


		double h,g,time1,f1,Df,dtime,time_0,xmid,xa,xb,count=0;
		double eps=std::numeric_limits<double>::epsilon(), tau=sqrt(eps); // square root of the machine epsilon is the tolerance for error
		 //we start Newton-Raphson method with finite difference
		while (fabs(Temp_fncb(r_0,r,r_1,fabs(velocity),Dm,beta,time,L,U,Temp_max))>tau){
			//cout << "After Temp_fncb" << endl;
			h=time*sqrt(eps);              //we start Newton-Raphson method with finite difference
			g=Temp_fncb(r_0,r,r_1,fabs(velocity),Dm,beta,time,L,U,Temp_max);
			time1=time+h;
			f1=Temp_fncb(r_0,r,r_1,fabs(velocity),Dm,beta,time1,L,U,Temp_max);
			Df=(f1-g)/h;                	//we use finite differenct to find Df in order to deal with the problem of diff(f)=0
			dtime = -g/Df;                 // calculating dtime,
			time_0=time;
			time = time + dtime;                 // updating the value of time
			count=count+1;
			xmid=(xmin+xmax)/2;         //here we have a bisection method to assist if the value of time_k+1 is out of bound under or over
			if (time<xmin || time>xmax){
				xa=xmid;
				if (Temp_fncb(r_0,r,r_1,fabs(velocity),Dm,beta,xmin,L,U,Temp_max)*Temp_fncb(r_0,r,r_1,fabs(velocity),Dm,beta,xmid,L,U,Temp_max)<0){
					xb=xmin;
				}
				else{
					xb=xmax;
				}
				xmin=min(xa,xb);
				xmax=max(xa,xb);
				time=(xmin+xmax)/2;
			}

			// modif VRG - 2017/30/03 -- changing the controls on the root finding algorithm to address particles being lost
			// These controls are new see ReturnDiffusionTime_final to see previous approaches.
			if (abs(time-time_0)/max(abs(time_0),1.)<tau && time<t_max){
				//cout << "Return Diffusion time - control 1 used " <<  endl;
				break;

			} else if(count>5000){
				//cout << "Return Diffusion time - control 2 used" <<  endl;
				time=ReturnDiffusionTime_bisection(r_0,r,r_1, fabs(velocity), beta,L,Dm,arr_time_max,U,Temp_max);
				break;
			}
		}

    double result=time;
	return result;

}

double FractureMesh::ReturnDiffusionTime_delta(vector<double> Time_Dist, vector<double> D_CDF, double U){
//double FractureMesh::ReturnDiffusionTime_delta(double *Time_Dist, vector<double> D_CDF, double U){
	double t_U;  int h4=D_CDF.size();
	//double max_D=std::max_element(D_CDF,D_CDF.size());
	int k=0;
	if (U >= D_CDF[h4-1]){
	    			t_U=Time_Dist[h4];
	    cout << "we are using an approximation in draw t delt because the CDF does not reach 1 - can l'hopital fix that" << endl;
	}else if (U<=D_CDF[0]) {
				t_U=Time_Dist[0];
	}else{
		//just write a for loop
		for(int i=1;i<h4;i++){
			if (D_CDF[i]>U ){
				k=i;
				i=h4;// break;
			}
		}
	       t_U=Time_Dist[k-1]+(Time_Dist[k]-Time_Dist[k-1])*((U-D_CDF[k-1])/(D_CDF[k]-D_CDF[k-1]));
	}
	   if (t_U<0){cout << "Warning in FractureMesh.c ReturnDiffusionTime_Con: diff time is negative"  << t_U << endl;   }
	return t_U;
}

// VRG march 2018

std::vector<double> FractureMesh::ReturnCDF_Delta(vector<double> PMF_Delta){
	vector<double> CDF_Delta;
	//cout << "End ReturnPMF_Delta" << endl;
	int h=PMF_Delta.size();
	CDF_Delta.reserve(h);
	CDF_Delta.assign(1,PMF_Delta[0]);
		for(int i=1;i<h;i++){
			CDF_Delta.push_back(PMF_Delta[i]+CDF_Delta[i-1]);
		}
		return CDF_Delta;
}

std::vector<double> FractureMesh::Shave_PMF(vector<double> PMF_k_end,int k){
	vector<double> sPMF;
	int h2= PMF_k_end.size();
	double PMF_k_end_1=PMF_k_end[k];
	sPMF.push_back(PMF_k_end_1);
	for(int i=k;i<h2-k+1;i++){
		sPMF.push_back(PMF_k_end[i+1]);
	}
	return sPMF;
}


// VRG November 2017  - Using Convolution and Generating function we generate a 3 distribution from two known distributions.
// We feed the function the Probability mass function of t_end and t_beg
//std::vector<double> FractureMesh::ReturnPMF_Delta(vector<double> PMF_end, vector<double> PMF_beg){
std::vector<double> FractureMesh::ReturnPMF_Delta(vector<double> PMF_End, vector<double> PMF_Beg){

    std::vector<double> PMF_Delta;
	double bD;
	int h=PMF_End.size();
	// prealocating memmory for PMF_delta
	PMF_Delta.reserve(h);
//  TRYING something new march 6 2018
/*				for (int k=0; k< h; k++){
					bD=0;
					if ( fabs (k-1)<1e-6){
						bD +=PMF_Beg[1]*PMF_Delta[k-1];
						cout << "bD"<< bD << "PMF_End[k]" <<PMF_End[k]<< "k="<< k<< "PMF_Beg[0]"<< PMF_Beg[0]<< endl;
					}
					for (int i=1; i<k; i++){
						bD +=PMF_Beg[i]*PMF_Delta[k-i];
					}
					double value1=(1/PMF_Beg[0])*(PMF_End[k] - bD);
					PMF_Delta.push_back(value1);
				}



*/
	//  TRYING something new march 7 2018
	double D_k;
	for (int k=0; k<h; k++){
		bD=0;
		for (int l=1; l<k+1;l++){
				bD+=PMF_Beg[l]*PMF_Delta[k-l];
		}
		D_k= (1/PMF_Beg[0])*(PMF_End[k]-bD);
		PMF_Delta.push_back(D_k);
	}

/* Commented out march 6
				for (int k=0; k< h; k++){
					bD=0;
					for (int i=1; i<k; i++){
					//for (int i=1; i<k+1; i++){
						bD +=PMF_Beg[i]*PMF_Delta[k-i];}


					double value1=(1/PMF_Beg[0])*(PMF_End[k] - bD);
					PMF_Delta.push_back(value1);
				}
*/

/*
				double S_PMF_Delta=0, pk=0;
				for (int i=0; i<PMF_Delta.size();i++){S_PMF_Delta=S_PMF_Delta+PMF_Delta[i];}
				if (abs(S_PMF_Delta-1)>10e-2){ cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): PMF_Delta does not sum to unity"<< endl;}
				for (int i=0; i<PMF_Delta.size();i++){
							if(PMF_Delta[i]<0){
									pk=pk+1;
									//cout << "i" << i<< endl;
							}
				}
				if(pk>0){

					cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): PMF_Delta has negative values"<< PMF_Beg[0]<< PMF_End[0] << endl;
				}

*/

		return PMF_Delta;



	/* 1 mars 2018
    std::vector<double> PMF_Delta;
	//std::vector<double> bD;
	double bD;
	int h=PMF_end.size();
	// prealocating memmory for PMF_delta
	PMF_Delta.reserve(h);
	//cout << "h " << h << endl;
	//bD.assign(1,0);



				//cout << "test5a h=" << h << endl;
				for (int k=0; k< h; k++){
					bD=0;
					for (int i=1; i<k+1; i++){
						bD +=PMF_beg[i]*PMF_Delta[k-i];
						//bD +=PMF_beg_test[i]*PMF_Delta_test[k-i];
						//bD+=PMF_beg[i];
						//bD+=1e-6;
					}

					double value1=(1/PMF_beg[0])*(PMF_end[k] - bD);
					PMF_Delta.push_back(value1);
					//PMF_Delta[k]=inv_pmf_beg_0*(PMF_end[k] - bD);
					//PMF_Delta_test[k]=inv_pmf_beg_0*(PMF_end_test[k] - bD)*0;
					//PMF_Delta[k]=bD;
					//else{PMF_Delta_test[k]=bD;cout << "test2"  << endl;}
				}
				//cout << "Er thetta sami kodi og i rvk 107 koppia yfir" << endl;
		return PMF_Delta;

	//clock_t t_begin=clock();
	 *
	 * */

/*
	int N=1e5;double value1;
	vector<double> vect1(N);
	for (int i=0;i<N;i++){
		value1=0.0;
		for (int j=1;j<i+1;j++){
			value1+=1e-6;
		}
		//cout << value1 << endl;
		//vect1[i]=value1;
	}

	//cout << "CPU Time = " << (clock()-t_begin)/CLOCKS_PER_SEC << endl;

	return vect1;
	*/
/*
	int N=PMF_end.size();
			std::vector<double> PMF_Delta(N);

			int h=PMF_end.size();

			float *PMF_Delta_test=(float*)calloc(PMF_end.size(),sizeof(float));
			float *PMF_end_test=(float*)calloc(PMF_end.size(),sizeof(float));
			float *PMF_beg_test=(float*)calloc(PMF_end.size(),sizeof(float));
			for (int k=0; k< h; k++){
				PMF_beg_test[k]=PMF_beg[k];
				PMF_end_test[k]=PMF_end[k];
			}

			double bD;
			// preallocating memory for PMF_delta
			//PMF_Delta.reserve(h);
			float inv_pmf_beg_0=1.0/PMF_beg[0];

			//cout << "test5a h=" << h << endl;
			for (int k=0; k< h; k++){
				bD=0.0;
				for (int i=1; i<k+1; i++){
					//bD +=PMF_beg[i]*PMF_Delta[k-i];
					//bD +=PMF_beg_test[i]*PMF_Delta_test[k-i];
					//bD+=PMF_beg[i];
					bD+=1e-6;
				}

				double value1=(1/PMF_beg[0])*(PMF_end[k] - bD);
				//PMF_Delta.push_back(value1);
				//PMF_Delta[k]=inv_pmf_beg_0*(PMF_end[k] - bD);
				//PMF_Delta_test[k]=inv_pmf_beg_0*(PMF_end_test[k] - bD)*0;
				PMF_Delta[k]=bD;
				//else{PMF_Delta_test[k]=bD;cout << "test2"  << endl;}
			}
			cout << "test5a after" << endl;
	return PMF_Delta;
*/
	/*  Thessi utgafa gefur retta nidurstodu en virkar ekki 1e6
	    std::vector<double> PMF_Delta;
		std::vector<double> bD;
		int h=PMF_end.size();
		// prealocating memmory for PMF_delta
		PMF_Delta.reserve(h);
		cout << "h " << h << endl;
		bD.assign(1,0);

		for (int k=0; k< h; k++){
		//	bD.assign(1,0);  	//	  THIS WORKS FOR 1E5 PROB 1E6  BUT it gives incorrect results so...
			double sum=0;   double geyma;

			for (int i=1; i<k+1; i++){

				double value2 =PMF_beg[i]*PMF_Delta[k-i];

				bD.push_back(value2);

				if (i>1){
					 geyma=bD[bD.size()-1];
					 sum+=geyma;
				}else{
					geyma=bD[bD.size()-1];
				    sum=geyma;
				}

			}
			double value1=(1/PMF_beg[0])*(PMF_end[k] - sum);
			PMF_Delta.push_back(value1);
		}
		cout << "h2 " << h << endl;
		return PMF_Delta;

/*
	/*
	std::vector<double> PMF_Delta;
	double bD;
	int h=PMF_end.size();
	// prealocating memmory for PMF_delta
	PMF_Delta.reserve(h);
	//cout << "h " << h << endl;
	//bD.assign(1,0);

	for (int k=0; k< h; k++){
		//bD.assign(1,0);
		bD=0;
		//double sum=0;   double geyma;

		for (int i=1; i<k+1; i++){

			//double value2 =PMF_beg[i]*PMF_Delta[k-i];
			//bD+=value2;
			//double value2 =
			bD+=PMF_beg[i]*PMF_Delta[k-i];


		}
		double value1=(1/PMF_beg[0])*(PMF_end[k] - bD);
		PMF_Delta.push_back(value1);
	}
	//cout << "h2 " << h << endl;
	return PMF_Delta;
	*/
}

std::vector<double> FractureMesh::ReturnCDF_Delta(double r_0,double r,double r_1,double Dm,double velocity,double advection_time, double t_min, double t_max,double sigma_end,double sigma_beg, double Nt_arr,vector<double>Time_Dist){
	std::vector<double> PMF_delta;
	std::vector<double> PMF_end,CDF_end;
	std::vector<double> PMF_beg,CDF_beg;
	// we need to check if the first value in the vector is to be shaved after the check the distributions are named:
	vector<double> PMF_Beg, PMF_End, PMF_Delta,CDF_Delta;
	int h= Time_Dist.size();

  	double Temp_max_end=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_max,sigma_end);
 	double Temp_max_beg=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_max,sigma_beg);
  	double advection_time_end=(sigma_end)/fabs(velocity);
  	double advection_time_beg=(sigma_beg)/fabs(velocity);
  	//std::cerr << "Before slow zone \n";
  	// G_X1pX2
	for (int i=0; i<h; i++){ CDF_end.push_back(Temp_fnc_c(r_0,r,r_1,Dm,fabs(velocity),advection_time_end,Time_Dist[i],sigma_end,Temp_max_end));}
	// G_X2
	for (int i=0; i<h; i++){ CDF_beg.push_back(Temp_fnc_c(r_0,r,r_1,Dm,fabs(velocity),advection_time_beg,Time_Dist[i],sigma_beg,Temp_max_beg));	}
	// From numerical CDF to numerical PMF
	PMF_beg.assign(1,CDF_beg[0]);
	PMF_end.assign(1,CDF_end[0]);
	for(int i=1;i<h;i++){
	    PMF_beg.push_back(CDF_beg[i]-CDF_beg[i-1]);
	    PMF_end.push_back(CDF_end[i]-CDF_end[i-1]);
	    //xg(i)=x(i);
	}


	/*
	 * int h2= PMF_beg.size();
	PMF_Beg[0]=PMF_beg[1]+PMF_beg[0];
	PMF_End[0]=PMF_end[1]+PMF_end[0];
	//double PMF_beg_1=PMF_beg[0]+PMF_beg[1];
	//PMF_Beg[0]=PMF_beg_1;
	//double PMF_end_1=PMF_end[0]+PMF_end[1];
	//PMF_End.assign(1,PMF_end_1);
	for(int i=1;i<h2-1;i++){
		//PMF_Beg[i]=(PMF_beg[i+1]);
		PMF_Beg.push_back(PMF_beg[i+1]);
		PMF_End.push_back(PMF_end[i+1]);
	  //   xg2=xg(2:end);
	}*/


/*
	double PMF_beg_1=PMF_beg[0]+PMF_beg[1];
	PMF_beg[1]=PMF_beg[1]+PMF_beg[0];
	PMF_end[1]=PMF_end[1]+PMF_end[0];
	for(int i=0;i<h2-1;i++){
		PMF_Beg.push_back(PMF_beg[i+1]);
		PMF_End.push_back(PMF_end[i+1]);
	  //   xg2=xg(2:end);
	}
	*/
	/*int kk=4;
	for (int i=0; i<kk;i++){
		PMF_beg_1=
	}*/
	/*
	for(int i=0; i<h2;i++){
		if (i<1){ PMF_Beg.push_back(PMF_beg_1);}

	}

		PMF_Beg[i]=PMF_beg[i]
*/
	//THIS IS NOT WORKING PMF_Beg and PMF_End do not match b_F and e_F in Matlab!!!!! lok at that
	int h2= PMF_beg.size();
	double PMF_beg_1=PMF_beg[0]+PMF_beg[1];
	double PMF_end_1=PMF_end[0]+PMF_end[1];
	//cout <<  PMF_Beg.size() << endl;
	PMF_Beg.push_back(PMF_beg_1);
	PMF_End.push_back(PMF_end_1);

	for(int i=1;i<h2-1;i++){
		PMF_Beg.push_back(PMF_beg[i+1]);
		PMF_End.push_back(PMF_end[i+1]);
	}
	/*
	int h9=PMF_Beg.size();
	cout <<"h2= PMF_beg.size();"<< h2<< "h9=PMF_Beg.size();" << h9 << endl;
	cout	<<"PMF_Beg[0]"<<  PMF_Beg[0] << "  PMF_Beg[1]"<<  PMF_Beg[1] << "PMF_Beg[2]"<<  PMF_Beg[2] << "PMF_Beg[h9-1]"<<  PMF_Beg[h9-1] << endl; //<<  "PMF_Beg[1]"<< PMF_Beg[1]<<  "PMF_Beg[2]"<< PMF_Beg[2]<< endl;
	cout	<<"PMF_Beg[3]"<<  PMF_Beg[3] << "  PMF_Beg[4]"<<  PMF_Beg[4] << "PMF_Beg[5]"<<  PMF_Beg[5] << "PMF_Beg[h9-1]"<<  PMF_Beg[h9-1] << endl; //<<  "PMF_Beg[1]"<< PMF_Beg[1]<<  "PMF_Beg[2]"<< PMF_Beg[2]<< endl;
	cout <<"PMF_end[0]"<<PMF_end[0] << "PMF_end[1]"<< PMF_end[1]<<"PMF_end[2]"<<  PMF_end[2]<<"PMF_end[h2-1]"<<  PMF_end[h2-1] << endl; //<<  "PMF_Beg[0]"<< PMF_Beg[0]<<"PMF_End[0]"<<  PMF_End[0]
*/

	//int rk0=PMF_Beg.size(); int rk2=PMF_End.size();
	//cout << "  rk0=" << rk0 <<"  PMF_Beg[rk-1]"<< PMF_Beg[rk0-1]<<"  rk2=" << rk2 <<"  PMF_End[rk2-1]"<<  PMF_End[rk2-1] << endl;

	/*
	int h2= PMF_beg.size();
	if (PMF_beg[0]<5){
		// we are shaving of the first vector values because of the l'hopital situation
		for(int i=1;i<h2;i++){
			PMF_Beg.push_back(PMF_beg[i]);
			PMF_End.push_back(PMF_end[i]);
		  //   xg2=xg(2:end);
		}
	 }else{
		 for(int i=0;i<h2;i++){
		 			PMF_Beg.push_back(PMF_beg[i]);
		 			PMF_End.push_back(PMF_end[i]);
		 		}
	 }
*/

	//  This is the original approach without the linear extension March 1 2018
		PMF_Delta=ReturnPMF_Delta(PMF_End, PMF_Beg);


		double S_PMF_Delta=0; int vk=0;
		for (int i=0; i<PMF_Delta.size();i++){S_PMF_Delta=S_PMF_Delta+PMF_Delta[i];}

		if (abs(S_PMF_Delta-1)>10e-2){ cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): PMF_Delta does not sum to unity"<< endl;}

			for (int i=0; i<PMF_Delta.size();i++){ if(PMF_Delta[i]<0){vk=vk+1;   } }
			//if(vk>30){cout << "tekk out the CDF the PMF has negatives" << endl;}
			/*		if(vk>0){

			cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): PMF_Delta has negative values"<< endl;
			int vk3=PMF_Delta.size();
			//cout << "PMF_Delta[0]"<< PMF_Delta[0] <<"PMF_Delta[1]"<< PMF_Delta[1] << "PMF_Delta[vk3-1]"<< PMF_Delta[vk3-1] <<"vk3="<< vk3-1 << endl;
			//cout << "CDF_beg[0]"<< CDF_beg[0] <<"CDF_end[0]"<< CDF_end[0] << "CDF_beg[1]"<< CDF_beg[1] <<"CDF_end[1]"<< CDF_end[1] <<"CDF_end[2]"<< CDF_end[2] << endl;
			int vk0=PMF_Beg.size(); int vk2=PMF_End.size();
			//cout << "  vk0=" << vk0 <<"  PMF_Beg[vk-1]"<< PMF_Beg[vk-1]<<"  vk2=" << vk2
			cout	<<"PMF_End[0]"<<  PMF_End[0] << "  PMF_End[1]"<<  PMF_End[1] << "PMF_End[2]"<<  PMF_End[2] << endl; //<<  "PMF_Beg[1]"<< PMF_Beg[1]<<  "PMF_Beg[2]"<< PMF_Beg[2]<< endl;
			cout <<"PMF_end[0]"<<PMF_end[0] << "PMF_end[1]"<< PMF_end[1]<<"PMF_end[2]"<<  PMF_end[2] << endl; //<<  "PMF_Beg[0]"<< PMF_Beg[0]<<"PMF_End[0]"<<  PMF_End[0]
			//PMF_beg[0]<< PMF_end[0] << "PMF_Delta[0] "<< PMF_Delta[0] <<" D"<< Dm << " u "<< fabs(velocity) <<" L_beg"<< sigma_beg <<" L_end"<< sigma_end <<"T_max_end "<< Temp_max_end  <<endl;
			//cout << "Time_Dist" << Time_Dist[0]<< Time_Dist[1] << "tmin" << t_min<< "tmax" << t_max<<"ntarr" << Nt_arr<< endl;
		}*/

/*
	//// 1 mars 2018 VRG: Linear extension
	// First for the action area we use the convolution as is for the first K elements of the vector.
	double K=3e5;
	std::vector<double>  PMF_End_K, PMF_Beg_K, PMF_Delta_K, PMF_Delta_S;

	for (int k=0; k< K; k++){ PMF_End_K.push_back(PMF_End[k]);}
	for (int k=0; k< K; k++){ PMF_Beg_K.push_back(PMF_Beg[k]);}

	PMF_Delta_K=ReturnPMF_Delta(PMF_End_K, PMF_Beg_K);

	// This is the actual linear extension
	int v=K-(Nt_arr-K)/K*0.3;
	double t_1=Time_Dist[v];
	double t_2=Time_Dist[K];
	double d_1=PMF_Delta_K[v];
	double d_2=PMF_Delta_K[K];
	double m= ((d_2-d_1)/(t_2-t_1));

	// defining PMF d_F_star
	//cout << "Before PMF_Delta_star" << endl;
	for (int i=0; i< K; i++){PMF_Delta.push_back(PMF_Delta_K[i]) ;}
	for (int i=K; i< Time_Dist.size()-1; i++){  PMF_Delta.push_back(d_2-m*(Time_Dist[i]-t_2));	}
	//cout << "After PMF_Delta_star " << endl;


*/
/*
	//cout << "End ReturnPMF_Delta" << endl;
	int h3=PMF_Delta.size();
	CDF_Delta.reserve(h);
	CDF_Delta.assign(1,PMF_Delta[0]);
		for(int i=1;i<h3;i++){
			CDF_Delta.push_back(PMF_Delta[i]+CDF_Delta[i-1]);
		}
		*/
		CDF_Delta= ReturnCDF_Delta(PMF_Delta);



		int pk2=0; int pop=0, pk21=0, pop21=0;
		for (int i=0; i<CDF_Delta.size();i++){
			if(CDF_Delta[i] <0){ pop=pop+1;}
			if(CDF_Delta[i] > 1.02){pk2=pk2+1;}
		}
		if(pk2>0){
						std::vector<double> PMF_end_s,PMF_beg_s, PMF_Delta_s, CDF_Delta_s;
						int s=2;
						PMF_end_s=Shave_PMF(PMF_end,s);
						PMF_beg_s=Shave_PMF(PMF_beg,s);
						PMF_Delta_s=ReturnPMF_Delta(PMF_end_s, PMF_beg_s);
						CDF_Delta_s=ReturnCDF_Delta(PMF_Delta_s);

						for (int i=0; i<CDF_Delta_s.size();i++){
							if(CDF_Delta_s[i] <0){ pop21=pop21+1;}
							if(CDF_Delta_s[i] > 1.02){pk21=pk21+1;}
						}
						if(pop21>0){cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): CDF_Delta_s has negative values"<< endl;
						}else if(pk21>0){cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): CDF_Delta_s has values greater then 1"<< endl;
						}else{
							CDF_Delta=CDF_Delta_s;
							//cout << "shave"<< endl;
						}
						/*
						cout <<"PMF_end_s[0]"<<PMF_end_s[0] << "PMF_end_s[1]"<< PMF_end_s[1]<<"PMF_end_s[2]"<<  PMF_end_s[2] << endl;
						cout <<"PMF_beg_s[0]"<<PMF_beg_s[0] << "PMF_beg_s[1]"<< PMF_beg_s[1]<<"PMF_beg_s[2]"<<  PMF_beg_s[2] << endl;
						cout <<"PMF_Delta_s[0]"<<PMF_Delta_s[0] << "PMF_Delta_s[1]"<< PMF_Delta_s[1]<<"PMF_Delta_s[2]"<<  PMF_Delta_s[2] << endl;
						cout <<"CDF_Delta_s[0]"<<CDF_Delta_s[0] << "CDF_Delta_s[1]"<< CDF_Delta_s[1]<<"CDF_Delta_s[2]"<<  CDF_Delta_s[2] << endl;
						int Yend=PMF_end_s.size(),Ybeg=PMF_beg_s.size(),YpDelta=PMF_Delta_s.size(),YcDelta=CDF_Delta_s.size();
						cout <<"PMF_end_s size "<<Yend<<"PMF_beg_size"<<Ybeg<<"PMF_Delta_s[YpDelta-1]"<<YpDelta<<"CDF_Delta_s[YcDelta-1]"<<YcDelta  << endl;
						cout <<"PMF_end_s size "<<PMF_end_s[Yend-1]<<"PMF_beg_size"<<PMF_beg_s[Ybeg-1]<<"PMF_Delta_s[YpDelta-1]"<<PMF_Delta_s[YpDelta-1]<<"CDF_Delta_s[YcDelta-1]"<<CDF_Delta_s[YcDelta-1]  << endl;
						cout <<"PMF_end_s[0]"<<PMF_end_s[0]<<"PMF_beg_s[0]"<<PMF_beg_s[0]<<"PMF_Delta_s[0]"<<PMF_Delta_s[0]<<"CDF_Delta_s[0]"<<CDF_Delta_s[0]  << endl;
						 */
	//					cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): CDF_Delta has values greater then 1"<< endl;
	//					cout << "velocity"<< velocity<< "sigma_end"<< sigma_end<<"sigma_beg"<<sigma_beg << endl;
						/*
						cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): CDF_Delta has values greater then 1"<< endl;
						cout << "CDF_beg[0]"<< CDF_beg[0] <<"CDF_end[0]"<< CDF_end[0] << "CDF_beg[1]"<< CDF_beg[1] <<"CDF_end[1]"<< CDF_end[1] << endl;
						cout <<" D"<< Dm << " u "<< fabs(velocity) <<" L_beg"<< sigma_beg <<" L_end"<< sigma_end <<"T_max_beg "<< Temp_max_beg <<"T_max_end "<< Temp_max_end  <<endl;
						cout << "Time_Dist" << Time_Dist[0]<< Time_Dist[1] << "tmin" << t_min<< "tmax" << t_max<<"Ntarr" << Nt_arr<< endl;
						int hk2= CDF_beg.size();
						int hk3= CDF_end.size();
						cout << "Time_Dist[1]" << Time_Dist[1] << "CDF_beg[0]"<< CDF_beg[0] <<"CDF_end[0]"<< CDF_end[0] <<endl;
						cout <<  "CDF_beg[1]"<< CDF_beg[1] <<"CDF_end[1]"<< CDF_end[1] <<endl;
						cout << "hk3" << hk3 << "hk2= CDF_beg.size()"<< hk2<< "CDF_beg[hk2-1] "<< CDF_beg[hk2-1] <<"CDF_end[hk3-1]"<< CDF_end[hk3-1] <<endl;


						 	cout <<  "CDF_Delta[0]" <<CDF_Delta[0] << "CDF_Delta[1]" <<CDF_Delta[1] <<"CDF_Delta[2]" <<CDF_Delta[2] <<"CDF_Delta[3]" <<CDF_Delta[3] <<endl;
						int hk=CDF_Delta.size();
						cout << "hk=CDF_Delta.size();" << hk << "CDF_Delta[hk-2]" <<CDF_Delta[hk-2] << "CDF_Delta[hk-1]" <<CDF_Delta[hk-1] <<endl;*/


		}else if (pop>0){
			cout << "WARNING in ReturnCDF_Delta (FractureMesh.cpp): CDF_Delta has negative values"<< CDF_Delta[0] << "CDF_Delt[1]" << CDF_Delta[1] << endl;
			cout << "pop" << pop  <<endl;
		}


	/*CDF_Delta.assign(1,0);
	for(int i=1;i<h3;i++){
		CDF_Delta.push_back(PMF_Delta[i-1]+CDF_Delta[i-1]);
	}*/

/*

	cout << "WRITTING CDF AND PMF IN FILES" << endl;
		string output_file1="./../Output/Results_pmf_beg.txt",
				output_file2="./../Output/Results_pmf_end.txt",
				output_file3="./../Output/Results_pmf_delta.txt",//,
				output_file4="./../Output/Results_cdf_delta.txt",
		        output_file5="./../Output/Results_Time_Dist.txt"; // vrg
		std::ofstream output1(output_file1.c_str(), std::ofstream::out),output2(output_file2.c_str(), std::ofstream::out), output3(output_file3.c_str(), std::ofstream::out), output4(output_file4.c_str(), std::ofstream::out), output5(output_file5.c_str(), std::ofstream::out);
			//,output4(output_file4.c_str(), std::ofstream::out);
		for (int i=0; i<PMF_Beg.size(); i++){
			output1 << PMF_Beg[i] << endl;
			output2 << PMF_End[i] << endl;
			output3 << PMF_Delta[i] << endl;
		}
		output1.close();output2.close();output3.close();

		for (int i=0; i<CDF_Delta.size(); i++){
			output4 << CDF_Delta[i] << endl;
		}

		for (int i=0; i<Time_Dist.size(); i++){
				output5 << Time_Dist[i] << endl;
			}
		output4.close();output5.close();
*/

	return CDF_Delta;
}

/*std::vector<double> FractureMesh::ReturnCDF_Delta(double r_0,double r,double r_1,double Dm,double velocity,double advection_time, double t_min, double t_max,double sigma_end,double sigma_beg, double Nt_arr,double *Time_Dist){
	cout << "begin ReturnCDF_Delta" << endl;
	std::vector<double> PMF_delta;
	std::vector<double> PMF_end,CDF_end;
	std::vector<double> PMF_beg,CDF_beg;
	// we need to check if the first value in the vector is to be shaved after the check the distributions are named:
	vector<double> PMF_Beg, PMF_End, PMF_Delta,CDF_Delta;
	//int h= Time_Dist.size();
	int h=Nt_arr;

  	double Temp_max_end=Temp_fnc_max(r_0,r,r_1,Dm,velocity,advection_time,t_max,sigma_end);
 	double Temp_max_beg=Temp_fnc_max(r_0,r,r_1,Dm,velocity,advection_time,t_max,sigma_beg);
  	double advection_time_end=(sigma_end)/fabs(velocity);
  	double advection_time_beg=(sigma_beg)/fabs(velocity);

  	// G_X1pX2
  	cout << "test1" << endl;
	for (int i=0; i<h; i++){ CDF_end.push_back(Temp_fnc_c(r_0,r,r_1,Dm,velocity,advection_time_end,Time_Dist[i],sigma_end,Temp_max_end));}
	cout << "test2" << endl;
	// G_X2
	for (int i=0; i<h; i++){ CDF_beg.push_back(Temp_fnc_c(r_0,r,r_1,Dm,velocity,advection_time_beg,Time_Dist[i],sigma_beg,Temp_max_beg));	}
	cout << "test3" << endl;
	// From numerical CDF to numerical PMF
	PMF_beg.assign(1,CDF_beg[0]);
	PMF_end.assign(1,CDF_end[0]);
	for(int i=1;i<h-1;i++){
	    PMF_beg.push_back(CDF_beg[i+1]-CDF_beg[i]);
	    PMF_end.push_back(CDF_end[i+1]-CDF_end[i]);
	    //xg(i)=x(i);
	}
	cout << "test4" << endl;
	int h2= PMF_beg.size();
	if (PMF_beg[0]<5){
		// we are shaving of the first vector values because of the l'hopital situation
		for(int i=1;i<h2;i++){
			PMF_Beg.push_back(PMF_beg[i]);
			PMF_End.push_back(PMF_end[i]);
		  //   xg2=xg(2:end);
		}
	 }else{
		 for(int i=0;i<h2;i++){
		 			PMF_Beg.push_back(PMF_beg[i]);
		 			PMF_End.push_back(PMF_end[i]);
		 		}
	 }
	cout << "test5" << endl;
	PMF_Delta=ReturnPMF_Delta(PMF_End, PMF_Beg);
	cout << "test6" << endl;
	int h3=PMF_Delta.size();
	CDF_Delta.reserve(h);
	CDF_Delta.assign(1,0);
	for(int i=1;i<h3;i++){
		CDF_Delta.push_back(PMF_Delta[i-1]+CDF_Delta[i-1]);
	}
	cout << "end ReturnCDF_Delta" << endl;
	return CDF_Delta;
}*/

// VRG Nov 2017  seg i > 1  function 3
double  FractureMesh::ReturnDiffusionTime_Con(double & res_time_seg,double advection_time,double Dm,double porosity,RngStream_a & rng_tracker,double t_in_fract,double sigma_beg,double previous_time, double Nt_arr,Declare_CDF_Map& Delta_CDF_Map,int pa_no){
	double t_max=MAX_TIME;// ;  // This is our observation time of the system i.e. 50 yearsm // OK
	double delta_time,diff_time;//,arr_time_max=t_max;
  	double U, beta, U_m, beta_m, r=(porosity*Dm)/(fabs(velocity)*aperture),r_0=r*r+1,r_1=pow(r,4.)+r*r+2;
     // Ensuring that min time requirement for the analytical solution holds
	double t_min=(10*10*10*10*aperture*aperture)/(porosity*porosity*Dm) *1.1;//4486;//
	std::vector<double> Delta_CDF;
	double sigma_end=sigma_beg+advection_time*fabs(velocity);
  	double Temp_max_end=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_max,sigma_end);
 	double Temp_max_beg=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_max,sigma_beg);
  	//Identifying the smallest temp value the analytical solution can handle based on its assumptions
  	double Temp_min_end=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_min,sigma_end)/Temp_max_end;
 	double Temp_min_beg=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_min,sigma_beg)/Temp_max_beg;
 	double Temp_min=max(Temp_min_end,Temp_min_beg);


    U = rng_tracker.uniform(0,1);
    if (U > Temp_min) {

    	   // Do a check to see if a Delta_CDF has been made for the same parameters (abs(sigma_end-sigma_beg), velocity)
    	   // If it has been then draw from there otherwise generate the Delta_CDF
    		double delta_t;	int h;
    	    	delta_t=(t_max-t_min)/(Nt_arr-1);
    	    std::vector<double> Time_Dist;
    		for (int i=0;i<Nt_arr;i++){Time_Dist.push_back(t_min+(i)*delta_t);}
    		h= Time_Dist.size();
    		/*double *Time_Dist=(double*)calloc(Nt_arr,sizeof(double));
    		for (int i=0;i<Nt_arr;i++){Time_Dist[i]=t_min+(i)*delta_t;}
    		h= Nt_arr;*/


int part_nr= pa_no;

    		//cout << "part_nr = " << part_nr +1<< endl;

    		if (!CDF_exist_and_return(fabs(velocity),fabs(sigma_end-sigma_beg),Delta_CDF_Map,Delta_CDF)){	// look for the CDF
    			// compute the cdf if it is not found in the map
    	    	  Delta_CDF=ReturnCDF_Delta( r_0,r,r_1,Dm,fabs(velocity),advection_time,t_min,t_max,sigma_end,sigma_beg,Nt_arr,Time_Dist);
    	    	  // store the computed cdf in the map
    	    	  // It has to be fabs(velcity, fabs(sigma_end),
    	    	  int particle_lim=7;
    	    	  if (pa_no<particle_lim){
    	    	  Delta_CDF_Map[make_pair(fabs(velocity),fabs(sigma_end-sigma_beg))]=Delta_CDF;
    	    	  	  	  if (pa_no >particle_lim-2){ cout << "We are one particle away from max storing and are at particle number "  << pa_no+1 << endl;	 }
    	    	  }

    		}

    		/// Time_Dist, Delta_CDF are correct so the question is what happens once it goes into the ReturnDiffusionTime_delta
    	  diff_time=ReturnDiffusionTime_delta(Time_Dist, Delta_CDF,U);


	  	   if (diff_time<0){ cout << "Warning in FractureMesh.c ReturnDiffusionTime_Con: diff time is negative"  << diff_time << endl;	  	   }



    }
    else{

		Get_Total_Time_From_HT1(U,res_time_seg,advection_time,Dm,porosity,previous_time);
		 diff_time=res_time_seg-advection_time;	// Non anymore we subtract advection_time
	}



    return diff_time;
}

// VRG November 2017 Using Convolution and Generating function  -- func 1
bool FractureMesh::Get_Total_Time_From_Advec_Time2D_Conv(double & res_time_seg,double advection_time,double Dm,double porosity,RngStream_a & rng_tracker,double t_in_fract,double sigma_beg, double previous_time, double Nt_arr,Declare_CDF_Map& Delta_CDF_Map,int pa_no){
	double t_max=MAX_TIME;// ;  // This is our observation time of the system i.e. 50 yearsm // OK
	double delta_time, res_diff_time_delta,res_diff_time_end, res_time_sigma;//,arr_time_max=t_max;
  	double sigma_end=sigma_beg+advection_time*fabs(velocity);

  	//cout << "sigma_beg at the beginning = " << sigma_beg << endl;

  	double U, beta, U_m, beta_m;
	double Dl=Dm; //9.16e-8;
	//if(fabs(Dl-Dm)>10e-9){ cout << "FractureMesh.C FYI D_l is not equal to D_t ensure that is the intention" << endl; }

	double r=(porosity*sqrt(Dm*Dl))/(fabs(velocity)*aperture),r_0=r*r+1,r_1=pow(r,4.)+r*r+2;
	// VRG March 10 2038 For the parallel validation figurei I had to have different longitudinal and transverse diffusion
	//double r=(porosity*Dm)/(fabs(velocity)*aperture),r_0=r*r+1,r_1=pow(r,4.)+r*r+2;


	// Modification Mai 23 2017 VRG for the 2D diffusion analytical solution to hold the following must be true:
	// Updated June 19 2017
	double R= (porosity*sqrt(Dm*Dl))/(fabs(velocity)*aperture);
	//double R= (porosity*sqrt(Dm*Dm))/(velocity*aperture);
	double K= ((porosity)/(2*aperture))*sqrt(Dm/Dl);
	double ter= (10*10*10*10*aperture*aperture)/(porosity*porosity*Dm);
	//cout << "R" <<  R<< endl;
	 if ( R <= 1){cout << "WARNING FOR analytical solution requires that R>1" << endl; cout << "velocity" <<  velocity << endl;exit(0);}
	 if ( K<= 1	  ){cout << "WARNING FOR analytical solution requires that K>1" << endl;}
	 if (	t_max <= ter	   ){cout << "WARNING FOR analytical solution requires that t_er>t_max" << endl;}

		 // July 11 2017 VRG ensuring that min time requirement for the analytical solution holds
	double t_min=(10*10*10*10*aperture*aperture)/(porosity*porosity*Dm) ;
	if (sigma_beg<EPSILON){
		//sigma_beg=EPSILON;
		//cout << " sigma_beg method H segm "<< sigma_beg << endl;
	  	double Temp_max_end=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_max,sigma_end);
	  	//cout << Temp_max_end<< endl;
	  	double Temp_min_end=Temp_fnc_max(r_0,r,r_1,Dm,fabs(velocity),advection_time,t_min,sigma_end)/Temp_max_end;
	  	double advection_time_end=(sigma_end)/fabs(velocity);
	  	//cout << "advection_time_end"<< advection_time_end << endl;
         U = rng_tracker.uniform(0,1);
         if (U<Temp_min_end){
        		//cout << "res_time_seg before" <<  res_time_seg << endl;
        	 	 Get_Total_Time_From_HT1(U,res_time_seg,advection_time,Dm,porosity,previous_time);
        	 	 // This creates an error I dunno why?
        	 	//res_time_seg=res_time_sigma;
         	 //	cout << "res_time_seg after" <<  res_time_seg << endl;
         }
         else{
        	 	 beta =  U* Temp_max_end;
        	 	 res_diff_time_end=ReturnDiffusionTime_HT2(r_0,r,r_1, fabs(velocity), beta,sigma_end,Dm,t_max,U,Temp_max_end);
        		//  	cout << "Temp_max_end" << Temp_max_end << endl;
        	 	//cout << "res diff time end" << res_diff_time_end << endl;
        	 	//cout << "advection time end" << advection_time_end << endl;
        	 	res_time_seg=res_diff_time_end+advection_time_end;
        	 	//cout << "res time seg " << res_time_seg << endl;
         }
		//	cout << "Test case for first segment" <<endl;
		//	res_time_seg=1;
     	//if (res_time_seg<0) {
     	//cout << "WARNING in (FractureMesh.cpp): diff_time is negative"<< sigma_beg  << endl;
     	//}
	}
	else{
		//	cout << " sigma_beg conv segm > "<< sigma_beg << endl;
		//	cou	t << "sigma_beg - is this a true messure if it is the first segment?" << sigma_beg<<endl; // I'm gettting 1e-10 for the firs seg
			double advection_time_seg=fabs(sigma_end-sigma_beg)/fabs(velocity);
			res_diff_time_delta=ReturnDiffusionTime_Con(res_time_seg,advection_time,Dm,porosity,rng_tracker,t_in_fract,sigma_beg,previous_time,Nt_arr,Delta_CDF_Map,pa_no);
		 //	cout << "res diff time delta" << res_diff_time_delta << endl;
		// 	cout << "advection time seg" << advection_time_seg << endl;
			res_time_seg = res_diff_time_delta+ advection_time_seg;
		 //	cout << "res time seg " << res_time_seg << endl;
	     	/*if (res_time_seg<0) {
	     	cout << "WARNING in (FractureMesh.cpp): diff_time is negative"<< sigma_beg  << "sigma_end" << sigma_end << "u" <<velocity << "t_in_frac" << t_in_fract << previous_time << endl;
	     	}*/
	}
	 //	cout << "res time seg " << res_time_seg << endl;

	//if (res_time_seg<0) {
//	cout << "WARNING in (FractureMesh.cpp): diff_time is negative"<< sigma_beg  << endl;
//	}
	return true;
}

bool FractureMesh::Get_Total_Time_From_HT1(double u,double & res_time,double advection_time,double Dm,double porosity,double previous_time){

//	cout << "U is" << u << endl;
//	cout << "res time in HT1 before " << res_time << endl;
//	cout << "previous_time" << previous_time << endl;

double B = sqrt(Dm)*porosity*advection_time/(this->aperture*0.5);
double value = boost::math::erfc_inv(u);
res_time = advection_time+0.25*B/value*B/value;
//cout << "res_time in HT1 after" << res_time << endl;
return true;

}

// Returns the diffusion time for transport process that assumes that we have 1D diffusion in the matrix
bool FractureMesh::Get_Total_Time_From_Advec_Time1D(double & res_time,double advection_time,double Dm,double porosity,RngStream_a & rng_tracker,double previous_time){

	/*
	cout << "phi" << porosity << endl;

	cout << "advection_timei" << advection_time << endl;
	cout << "Dm" << Dm << endl;
	cout << "aperture" << this->aperture << endl; */

	double B = sqrt(Dm)*porosity*advection_time/(this->aperture*0.5);
	double u = rng_tracker.uniform(0,1);
	double value = boost::math::erfc_inv(u);
	res_time = advection_time+0.25*B/value*B/value;
	return true;

}

/*bool FractureMesh::Get_Total_Time_From_Advec_Time(double & res_time,double advection_time,double Dm,double porosity,RngStream_a & rng_tracker,double t_in_fract,double L_in_fract,double Nt_arr,Declare_CDF_Map& Delta_CDF_Map,int pa_no){

// VRG Nov 2017 Convolution with Generating functions
	return Get_Total_Time_From_Advec_Time2D_Conv(res_time,advection_time,Dm,porosity,rng_tracker,t_in_fract,L_in_fract,t_in_fract, Nt_arr, Delta_CDF_Map,pa_no);

	// Produce the time if there is 1D diffusion in the matrix
	//return Get_Total_Time_From_Advec_Time1D(res_time,advection_time,Dm,porosity,rng_tracker,t_in_fract);
}*/

/** returns the advection time required for the scheme stability (transfer probability lower than transfer_proba_max)*/
double FractureMesh::Advection_Time_Computation(pointcpp<double> current_position, const Projection_Map & trans_positions, const double porosity, const double diffusion, const double transfer_proba_max, const double t_advec_end){

	//case without restriction
	if (transfer_proba_max == -1) return t_advec_end;
	//0. PARAMETERS
	if ((trans_positions.size()<0)||(trans_positions.size()>2)){
		std::cout<<"\n PB IN Advection_Time_Computation (Utilitary_Transfer.cpp) : res not defined";
		return -1;
	}

	//1. NO SURROUNDING FRACTURES
	if (trans_positions.size()==0) return t_advec_end;

	//2. ONE OR TWO SURROUNDING FRACTURES
	//barriers parameters
	double distance1 = 0, distance2 = 0;
	if (trans_positions.size()==1){
		distance1 = trans_positions.begin()->first.Point_Distance(current_position);
		distance2 = distance1;
	}
	else if (trans_positions.size()==2){
		Projection_Map::const_iterator it = trans_positions.begin();
		distance1 = it->first.Point_Distance(current_position);
		it++;
		distance2 = it->first.Point_Distance(current_position);
	}
	double mean_transfer_time = distance1*distance2/(2*diffusion);
	double advection_time = aperture*std::sqrt(mean_transfer_time)*boost::math::erfc_inv(1-transfer_proba_max)/(porosity*std::sqrt(diffusion));
	return std::min(advection_time, t_advec_end);
}

pointcpp<double> FractureMesh::Mesh_Scale_Advection(pointcpp<double> M_in, double dt, bool opposite){
	pointcpp<double> v_inter = velocity_interpolation(M_in);
	if (opposite) v_inter = pointcpp<double>(-v_inter.i,-v_inter.j,-v_inter.k);
	return M_in + v_inter * dt;
}

/// Gets the interpolated velocity at point M
pointcpp<double> FractureMesh::velocity_interpolation(const pointcpp<double> & M){
	pointcpp<double> res(0.,0.,0.);
	CgalVector2D vect(p_ori.p, p_tar.p);
	res.i = vect.x()/CGAL::sqrt(vect.squared_length()) * velocity;
	res.j = vect.y()/CGAL::sqrt(vect.squared_length()) * velocity;
	return res;
}

int FractureMesh::Is_Collinear(CgalPoint2D pt){
	if (CGAL::collinear(p_ori.p,p_tar.p,pt)) return 1;
	return 0;
}

//returns 1 if the current mesh intersects the segment seg
//inter is the intersection point
int FractureMesh::Is_Intersected(Segment2D & seg, CgalPoint2D & inter){
	Segment2D seg_mesh(p_ori.p,p_tar.p);
	if (seg_mesh.intersection(seg,inter)) return 1;
	return 0;
}

Segment2D FractureMesh::ReturnSegment(){
	Segment2D seg_mesh(p_ori.p,p_tar.p);
	return seg_mesh;
}

// fracture length
double FractureMesh::ReturnLength(){
	return CGAL::sqrt(CGAL::squared_distance(this->p_ori.p,this->p_tar.p));
}

double FractureMesh::ReturnConductance(){
	return RHO*G*std::pow(aperture,2)/(12*MU);
}
