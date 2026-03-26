/*
 * Transfer.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */


#include "Transport.h"
#include "../Utilitaries/Laplace_Stehfest.h"
#include "../Utilitaries/scale.h"
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;


/**
* transfer = true if the particle is transfered from current_position to one of the surrounding barriers before the time t_max
* trans_positions : possible positions of transfer
* returns the mesh where the particle must be transfered (next_mesh)
* returns the required time to transfer (transfer_time)
*/
bool Transport::Transfer_Probability_and_Transfer_Time(Projection_Map trans_positions, const double t_max, int & next_mesh, double & transfer_time,pointcpp<double> & M_proj){
	//0. PARAMETERS
	//0.1. outputs
	next_mesh = -1; transfer_time = -1;
	if (trans_positions.size()==0){return false;}
	//0.2. barrier definition
	double coeff = 1;
	std::map<int,double> barrier_positions; std::pair<double,double> distances = std::make_pair(-1,-1);
	for (Projection_Map::const_iterator it = trans_positions.begin(); it != trans_positions.end(); it++){
		double dist = coeff*it->second.second;
		barrier_positions[it->second.first] = dist;
		if (distances.first == -1) distances.first = dist;
		else distances.second = dist;
	}
	//1. TRANSFER PROBABILITY
	//1.0. determination of the studied barrier
	double u_simu = rng_tracker.uniform(0,1);
	//1.1. computation of the maximal transfer probability
	double u_max = 1;
	if (!Transfer_Probability_Computation(distances, t_max, u_max)){
		std::cout<<"\n proba_transfer not defined: u_max = " << u_max << endl;
		cout << distances.first << " " << distances.second << " " << t_max << endl;
	}
	//1.2. determination of the transfer
	if (u_simu>u_max){return false;}//no transfer
	else{	//transfer time determination
		if (!Transfer_Time_Computation(distances, u_simu, transfer_time)){std::cout<<"\n transfer_time not defined";}
		//next_mesh determination
		u_simu = u_simu/u_max;
		//if (!Barrier_Choice(barrier_positions, u_simu, next_mesh)){std::cout<<"\n PB IN Transfer_Probability_and_Transfer_Time (Tracker_Transfer.cpp) : mesh probabilities not defined";}
		std::string next_position;
		if (!Barrier_Choice(distances, u_simu, next_position)){std::cout<<"\n PB IN Transfer_Probability_and_Transfer_Time (Transfer.cpp) : mesh probabilities not defined";}
		//transfer_time = 0.5*transfer_time;	// Sudicky hypothesis
		// next position
		Projection_Map::iterator it_choice=trans_positions.begin();
		if (next_position=="first"){M_proj=it_choice->first;next_mesh=it_choice->second.first;}
		else if (next_position=="second"){it_choice++;M_proj=it_choice->first;next_mesh=it_choice->second.first;}
		else {std::cout<<"\n PB IN Transfer_Probability_and_Transfer_Time (Transfer.cpp) : next position not defined";}
		return true;
	}
	return false;
}

/**
* returns the next mesh (mesh of transfer) from the transfer probability u
* returns distances(distance1,distance2) with distance1 the reached barrier
* u_simu modfication to take into account the arrival barrier
*/
bool Barrier_Choice(const std::map<int,double> & barrier_positions, const double u_simu, int & next_mesh){
	//0. PARAMETERS
	next_mesh = -1;
	//1. ONE NEIGHBOURING FRACTURE
	if (barrier_positions.size()==1){next_mesh = barrier_positions.begin()->first;}
	//2. TWO NEIGHBOURING FRACTURES
	else if (barrier_positions.size()==2){
		//2.1. distance definition
		std::map<int,double>::const_iterator it1 = barrier_positions.begin(), it2 = it1; it2++;
		double distance1 = std::fabs(it1->second);
		double distance2 = std::fabs(it2->second);
		//2.2. probability affectation
		double proba1 = distance2/(distance1+distance2);
		//probability of reaching the barrier at the distance 'distance1'
		if (u_simu<proba1) next_mesh = it1->first;
		else next_mesh = it2->first;
	}
	// 3. NOT IMPLEMENTING CASES
	else{
		cout<<"\n PB IN Barrier_Choice (Transfer.cpp) : case not implemented";
		cout << "Number of barriers = " << barrier_positions.size() << endl;
		return false;
	}
	return true;
}

bool Barrier_Choice(std::pair<double,double> distances, const double u_simu, std::string & next_position){
	//0. PARAMETERS
	next_position = NO_POSITION;
	//1. ONE NEIGHBOURING FRACTURE
	if (distances.first!=-1&&distances.second==-1){next_position = "first";}
	//2. TWO NEIGHBOURING FRACTURES
	else if (distances.first!=-1&&distances.second!=-1){
		//2.1. distance definition
		double distance1 = distances.first;
		double distance2 = distances.second;
		//2.2. probability affectation
		double proba1 = distance2/(distance1+distance2);
		//probability of reaching the barrier at the distance 'distance1'
		if (u_simu<proba1){next_position = "first";}
		else{next_position = "second";}
	}
	// 3. NOT IMPLEMENTING CASES
	else{
		cout<<"\n PB IN Barrier_Choice (Transfer.cpp) : case not implemented with distances";
		return false;
	}
	return true;
}

/**
* transfer_proba = the probability of transfer from current_position to one of the trans_positions at the time t_max
*/
bool Transport::Transfer_Probability_Computation(const std::pair<double,double> & distances, const double t_max, double & transfer_proba){
	//0. PARAMETERS
	double diffusion = phys_param.Dm;
	//1. NO POSSIBILITY OF TRANSFER
	if ((distances.first == -1)&&(distances.second == -1)) transfer_proba = 0;
	//2. ONE POSSIBILITY OF TRANSFER ->	FIRST PASSAGE TIME DISTRIBUTION
	else if (distances.second==-1){
		transfer_proba = Transfer_Probability_FPTD(distances.first, diffusion, t_max);
	}
	//3. TWO POSSIBILITIES OF TRANSFER -> TIME TRANSFER DISTRIBUTION
	else if ((distances.first!=-1)&&(distances.second!=-1)){
		return Transfer_Probability_Feller(distances, diffusion, t_max, transfer_proba);
	}
	else{
		std::cout<<"\n PB in Transfer_Probability_Computation (Transfer.cpp) : case not implemented";
		return false;
	}
	return true;
}

/**
* returns time required to transfer with the transfer probability u_simu
*/
bool Transport::Transfer_Time_Computation(const std::pair<double,double> & distances, const double u_simu, double & transfer_time){
	//0. PARAMETERS
	transfer_time = -1;
	double diffusion = this->phys_param.Dm;

	//1. ONE POSSIBILITY OF TRANSFER ->	FIRST PASSAGE TIME DISTRIBUTION
	if (distances.second==-1){
		transfer_time = Transfer_Time_FPTD(distances.first, diffusion, u_simu);
		return true;
	}
	//2. TWO POSSIBILITIES OF TRANSFER -> TIME TRANSFER DISTRIBUTION
	else if ((distances.first!=-1)&&(distances.second!=-1)){
		return Transfer_Time_Feller(distances, diffusion, u_simu, rng_tracker, transfer_time, Transfer_Time_Distribution);
	}

	std::cout<<"\n PB in Transfer_Probability (Transfer.cpp) : case not implemented";
	return false;
}

/**returns time required to reach the distance with the transfer probability u_simu*/
double Transfer_Time_FPTD(const double distance, const double diffusion, const double u_simu){
	double A = boost::math::erf_inv(1-u_simu);
	return std::pow(distance/A,2.)/(4*diffusion);
}

/**returns transfer_time : required time to reach one of the two barriers assuming that the transfer is done with the probability u_simu*/
bool Transfer_Time_Feller(const std::pair<double,double> & distances, const double diffusion, const double u_simu, RngStream_a & rng, double & transfer_time, Distribution & Transfer_Time_Distribution){

	//1. TIME COMPUTATION
	//1.1. time distribution determination
	//parameters
	std::pair<std::vector<double>,std::vector<double> > Single_Time_Distribution;

	//look for the required distribution in Transfer_Time_Distribution
	if ((!exists_key(Transfer_Time_Distribution,distances,Single_Time_Distribution,EPSILON_DIST,true))||
		(Single_Time_Distribution.second[0]>u_simu)||(Single_Time_Distribution.second[Single_Time_Distribution.second.size()-1]<u_simu)){
			if (!Transfer_Time_Distribution_Computation(Transfer_Time_Distribution, distances, diffusion, u_simu, rng)){
				std::cout<<"\n PB IN Transfer_Time_Feller (Transfer.cpp) : Single_Time_Distribution not defined";
				return false;
			}
	
			if (!exists_key(Transfer_Time_Distribution,distances,Single_Time_Distribution,EPSILON_DIST,true)){
				std::cout<<"\n PB IN Transfer_Time_Feller (Transfer.cpp) : Single_Time_Distribution not found";
				return false;
			}
	}

	//1.2. transfer time determination
	return Get_Time_Feller(Single_Time_Distribution, u_simu, transfer_time);
}

/**
* determines the required time to reach one of the barriers with the transfer probability u_simu
* Single_Time_Distribution : pair<time,transfer probability>
*/
bool Get_Time_Feller(std::pair<std::vector<double>,std::vector<double> > & Single_Time_Distribution, const double u_simu, double & transfer_time){

	//1. TRANSFER TIME DISTRIBUTION
	std::vector<double> t_real = Single_Time_Distribution.first;
	std::vector<double> transfer_proba = Single_Time_Distribution.second;

	//2. DETERMINATION OF THE CORRESPONDING TIME
	//2.1. look for the fixed probability
	int u_index = -1;
	if (!find_bound_in_vector(transfer_proba, u_simu, u_index)){
		std::cout << "\n PB IN Feller_Transfer_Time_Distribution (Transfer.cpp) : u not found" << endl;
                cout << "u = " << u_simu << endl;
               print(transfer_proba);
		 exit(0);
		return false;
	}
	//2.2. interpolation of the corresponding time
	transfer_time = linear_interpolation_2points(transfer_proba[u_index], transfer_proba[u_index+1], t_real[u_index], t_real[u_index+1], u_simu);

	return true;
}

/**
* returns false if the computation fails
* computation of the transfer time required to reach the distance1 before the distance2 given the transfer probability u_simu
*/
bool Transfer_Time_Distribution_Computation(Distribution & Transfer_Time_Distribution, const std::pair<double,double> & distances, const double diffusion, const double u, RngStream_a & rng){

	//0. PARAMETERS
	//check the current Transfer_Time_Distribution
	std::pair<std::vector<double>,std::vector<double> > Single_Time_Distribution;
	if ((exists_key(Transfer_Time_Distribution,distances,Single_Time_Distribution,EPSILON_DIST,true))&&
		(Single_Time_Distribution.second[0]>=u)&&
		(Single_Time_Distribution.second[Single_Time_Distribution.second.size()-1]<=u)) return true;

	//1. TIME MAX DETERMINATION
	double t_max = -1;
	if (!Feller_Maximal_Time_Determination(t_max, u, distances, diffusion)){
		std::cout<<"\n PB IN Transfer_Time_Distribution_Computation (Transfer.cpp) : t_max not defined";
		return false;
	}

	//2. DISTRIBUTION DETERMINATION
	//time discretization
	int Nt = 100;
	scale time(0,t_max,Nt+1,1);
	vector<double> t_real = time.scale_vector();
	t_real.erase(t_real.begin());
	std::vector<double> transfer_proba(t_real.size());
	//parameters for Stehfest inversion
	int Np = 8;	// Number of sample in the interval in Laplace space
	std::vector<double> Vi = Stehfest_Calcul_Vi(Np);

	//loop on time for computation
	for (size_t i = 0; i<t_real.size(); i++){
		// 1- SAMPLING IN LAPLACE SPACE
		std::vector<double> lambda = Stehfest_Sampling(Np, t_real[i]), p_laplace(Np+1);
		// 2- FUNCTION DETERMINATION IN LAPLACE SPACE
		double y = 0;
		for( int k=1; k<=Np; k++){
			if (!Feller_Analytical_Solution(y, lambda[k], distances, diffusion)){
				std::cout<<"\n PB IN Transfer_Time_Distribution_Computation (Transfer.cpp) : analytical solution not defined";
				return false;
			}
			p_laplace[k] = y;
		}
		// 3- INVERSION
		transfer_proba[i] = Stehfet_Inverse_Laplace(t_real[i], p_laplace, Np, Vi);
	}
	//addition of the case time=0 and correction
	for (size_t i = 0; i<transfer_proba.size(); i++)
		if (transfer_proba[i]<0) transfer_proba[i] = 0;
	t_real.insert(t_real.begin(),0); transfer_proba.insert(transfer_proba.begin(),0);

	//3. AFFECTATION
	/*//remove the existing value (if the distribution exists for a lower maximum time)
	remove_from_key(Transfer_Time_Distribution,distances,EPSILON_DIST,true);
	//add or replace the computed distribution
	Transfer_Time_Distribution[distances] = std::make_pair<std::vector<double>,std::vector<double> >(t_real,transfer_proba);*/

	// replace the existing distribution with the new one
	if (!replace_distribution(Transfer_Time_Distribution,distances,EPSILON_DIST,true,make_pair(t_real,transfer_proba))){
		// if there is no existing dist - create a new one
		//Transfer_Time_Distribution[distances] = std::make_pair<std::vector<double>,std::vector<double> >(t_real,transfer_proba);
		Transfer_Time_Distribution[distances] = std::make_pair(t_real,transfer_proba);
	}
	return true;
}

/**determines the maximal time of the Feller distribution given the probability of transfering to distance1 and distance2*/
bool Feller_Maximal_Time_Determination(double & t_max, const double u_simu, const std::pair<double,double> & distances, const double diffusion){
	//parameters
	int Np = 8;	// Number of sample in the interval in Laplace space
	std::vector<double> Vi = Stehfest_Calcul_Vi(Np);
	double time = EPSILON; double transfer_proba = -1;

	//look for the time corresponding to transfer_proba=u_max
	while (transfer_proba<u_simu){
		time = 10*time;
		// 1- SAMPLING IN LAPLACE SPACE
		std::vector<double> lambda = Stehfest_Sampling(Np, time), p_laplace(Np+1);
		// 2- FUNCTION DETERMINATION IN LAPLACE SPACE
		double y = 0;
		for( int k=1; k<=Np; k++){
			if (!Feller_Analytical_Solution(y, lambda[k], distances, diffusion)){
				std::cout<<"\n PB IN Feller_Transfer_Time_Distribution (Transfer.cpp) : analytical solution not defined";
				return false;
			}
			p_laplace[k] = y;
		}
		// 3- INVERSION
		transfer_proba = Stehfet_Inverse_Laplace(time, p_laplace, Np, Vi);
	}
	t_max = time;
	return true;
}


/**returns the probability of transfer to distance before the time T*/
double Transfer_Probability_FPTD(const double distance, const double diffusion, const double T){
	return 1-boost::math::erf(distance/(2*std::sqrt(diffusion*T)));
}

/**
* distance = pair(distance of barrier1, distance of barrier2)
* one of the two distance must be negative and the other one positive
* transfer_proba = the probability to reach one of the two barriers at the time T
*/
bool Transfer_Probability_Feller(const std::pair<double,double> & distances, const double diffusion, const double T, double & transfer_proba){

	//1. INVERSION OF THE PROBABILITY TO REACH THE FIRST BARRIER WITHOUT PASSING BY THE SECOND ONE

	// 1.0- Parameters
	//parameters for Stehfest inversion
	int Np = 8;	// Number of sample in the interval in Laplace space
	std::vector<double> Vi = Stehfest_Calcul_Vi(Np);

	// 1.1- Sampling parameter vector and extraction of parameter "g" according to Sudicky
	std::vector<double> lambda = Stehfest_Sampling(Np, T), p_laplace(Np+1);

	// 1.2- Sampling of function in Laplace space
	double y = 0;
	for( int k=1; k<=Np; k++){
		if (!Feller_Analytical_Solution(y, lambda[k], distances, diffusion)){
			std::cout<<"\n PB IN Feller_Transfer_Time_Distribution (Transfer.cpp) : analytical solution not defined";
			return false;
		}
		p_laplace[k] = y;
	}

	// 1.3- inversion and correction
	transfer_proba = Stehfet_Inverse_Laplace(T, p_laplace, Np, Vi);
	if (transfer_proba<0) transfer_proba = 0;

	return true;
}

/**
* probability to reach the position distance1 and the position distance2
* with the diffusion coefficient D in the laplace space
* the two distances must be positive
*/
bool Feller_Analytical_Solution(double & y, const double x, const std::pair<double,double> & distances, const double D){

	//PARAMETERS
	double distance1 = distances.first, distance2 = distances.second;

	//COMPUTATION
	double y1 = 0, y2 = 0;
	if ((!Feller_Analytical_Solution_Single(y1, x, distance1, distance2, D))||
		(!Feller_Analytical_Solution_Single(y2, x, distance2, distance1, D))){
		std::cout<<"\n PB IN Feller_Analytical_Solution (Transfer.cpp) : single analytical solution not defined";
		return false;
	}
	y = y1+y2;
	return true;
}

/**
* y = Feller analytical solution for the variable x
* probability to reach the barrier at the distance distance1 without passing by the barrier at the position distance2
* with the diffusion coefficient D in the laplace space
* the two distances must be positive
* it assumes that the lower barrier is at the distance distance1 and the upper one is at the distance distance2
*/
bool Feller_Analytical_Solution_Single(double & y, const double x, const double distance1, const double distance2, const double D){
	//barrier positions
	double d1 = -distance1, d2 = distance2;	//lower and upper position
	if (d1*d2>=0){
		std::cout<<"\n PB IN Feller_Analytical_Solution_Single (Transfer.cpp) : barrier positions not defined";
		return false;
	}

	//solution computation
	y = std::exp(std::sqrt(x/D)*d1)*(1-std::exp(-2*d2*std::sqrt(x/D)))/(x*(1-std::exp(2*(d1-d2)*sqrt(x/D))));

	//test of the solution
	if ((boost::math::isnan(y))||(!boost::math::isfinite(y))){
		std::cout<<"\n PB IN Feller_Analytical_Solution_Single (Transfer.cpp) : res not defined" << endl;
		return false;
	}
	return true;
}


