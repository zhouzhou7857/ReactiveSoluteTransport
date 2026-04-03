/*
 * Structures.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: delphine
 */


#include "Structures.h"
#include <iostream>
#include <utility>

using namespace std;


double get_min_value(std::map<int,double> values){
	double min_value = values.begin()->second;
	for (std::map<int,double>::iterator it=values.begin();it!=values.end();it++){
		if (it->second<min_value) min_value = it->second;
	}
	return min_value;
}

double get_min_value_pos(std::map<int,double> values){
	double min_value=1e20;
	for (std::map<int,double>::iterator it=values.begin();it!=values.end();it++){
		if (it->second>=0&&it->second<min_value){min_value = it->second;}
	}
	return min_value;
}

double get_max_value(std::map<int,double> values){
	double max_value = values.begin()->second;
	for (std::map<int,double>::iterator it=values.begin();it!=values.end();it++){
		if (it->second>max_value) max_value = it->second;
	}
	return max_value;
}

void print(std::map<int,double> values){
	for (std::map<int,double>::iterator it=values.begin(); it!=values.end(); it++){
		cout << it->first << " " << it->second << endl;
	}
}

void print(std::vector<double> values){
	int N = values.size();
	for (int i=0;i<N;i++){
		cout << values[i] << endl;
	}
}

void print(Distribution values){
	for (Distribution::iterator it=values.begin();it!=values.end();it++){
		cout << "(" << it->first.first << "," << it->first.second << "):" << endl;
		cout << "vector1:" << endl;
		int size1 = it->second.first.size();
		for (int i=0;i<size1;i++){cout << it->second.first[i] << endl;}
		cout << "vector2:" << endl;
		int size2 = it->second.second.size();
		for (int i=0;i<size2;i++){cout << it->second.second[i] << endl;}
	}
}

double define_delta(double min, double max, int N, int lin_log){
	if (N==1)
		return 0.;
	else if (lin_log==SCALE_LOG){
		if (min==0) min = 1;
		return pow(max/min,1./(N-1));
	}
	else if (lin_log==SCALE_LIN){
		return (max-min)/(N-1);
	}
	else{cout << "WARNING in define_delta (Results.cpp): scale option not defined" << endl;}
	return 0;
}

vector<double> define_scale(double min, double max, int N, int lin_log){

	vector<double> scale(N);
	double delta = define_delta(min,max,N,lin_log);
	for (int i=0;i<N;i++){
		if (lin_log==SCALE_LIN){
			scale[i] = min+delta*i;
		}
		else if (lin_log==SCALE_LOG){
			scale[i] = min*pow(delta,i);
		}
		else{cout << "WARNING in define_scale (Results.cpp): scale option not defined" << endl;}
	}
	return scale;
}

bool find_bound_in_vector(const std::vector<double> & vect, const double value, int & index){
	index = -1;
	for (size_t i = 0; i<vect.size()-1; i++){
		if ((vect[i]<=value)&&(vect[i+1]>=value)){
			index = i;
			return true;
		}
	}
	std::cout<<"\n PB IN find_bound_in_vector (Vector_Utilitary.cpp) : value not found";
	return  false;
}

double linear_interpolation_2points(const double xmin, const double xmax, const double ymin, const double ymax, const double x){
	return ymin + (x - xmin) / (xmax - xmin) * (ymax - ymin);
}

/**
* return true if value is found as the first value of one of the pairs in projection_map
* return the corresponding point
*/
bool find_projected_value(map<pointcpp<double>,pair<int,double> > proj_map, int value, pointcpp<double> & point){
	for (map<pointcpp<double>,pair<int,double> >::iterator it = proj_map.begin(); it!=proj_map.end(); it++){
		if (it->second.first == value){
			point = it->first;
			return true;
		}
	}
	cout << "WARNING in find_projected_value (Structures.cpp): projected value not found" << endl;
	return false;
}

pair<int,int> return_indices(int index,int Ny){return std::make_pair(floor(index/Ny),index-floor(index/Ny)*Ny);}

void print(Projection_Map projections){
	cout << "Projections: " << endl;
	for (Projection_Map::iterator it=projections.begin();it!=projections.end();it++){
		cout << "Point (" << it->first.i << "," <<  it->first.j << ") on mesh " << it->second.first << " at the distance " << it->second.second << endl;
	}
}

bool replace_distribution(std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double> > > & M, const std::pair<double,double> & key,const double epsilon,bool reverse,std::pair<std::vector<double>,std::vector<double> > dist){
	std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double> > >::const_iterator it;
	for(it=M.begin(); it!=M.end(); it++){
		if ((std::fabs(it->first.first-key.first)<epsilon)&&(std::fabs(it->first.second-key.second)<epsilon)){
			M[std::make_pair(it->first.first,it->first.second)]=dist;
			return true;
		}
	}
	if (reverse){
		if (replace_distribution(M, std::make_pair(key.second,key.first),epsilon,false,dist)){return true;}
	}
	return false;
}

bool exists_and_replace(MapSegmentTimesCum& map_seg_times,pair<pointcpp<double>,pointcpp<double> > pair_ext,map<double,double> map_times){
	bool find_key=false;map<double,double> current_map_times;
	for (MapSegmentTimesCum::iterator it1=map_seg_times.begin();it1!=map_seg_times.end();it1++){
		if (pair_ext.first.identic(it1->first.first,EPSILON)&&pair_ext.second.identic(it1->first.second,EPSILON)){
			find_key=true;current_map_times=it1->second;
			for (map<double,double>::iterator it2=map_times.begin();it2!=map_times.end();it2++){
				current_map_times[it2->first]=current_map_times[it2->first]+it2->second;
			}
			map_seg_times[make_pair(it1->first.first,it1->first.second)]=current_map_times;
			break;
		}
	}
	return find_key;
}

bool CDF_exist_and_return(double velocity,double diff_sigma,Declare_CDF_Map Delta_CDF_Map,vector<double>& Delta_CDF){
	for (Declare_CDF_Map::iterator it=Delta_CDF_Map.begin();it!=Delta_CDF_Map.end();it++){
		if ((fabs(it->first.first-velocity)<=EPSILON)&&(fabs(it->first.second-diff_sigma)<=EPSILON)){
			Delta_CDF=it->second;
			return true;
		}
	}
	return false;
}

