/*
 * Results.cpp
 *
 *  Created on: Jun 29, 2012
 *     Author: delphine
 */


#include "Results.h"
#include "../Utilitaries/Structures.h"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;



// determine the cumulative distribution of arrival times
Results::Results(map<int,double> arrival_times_,Parameters param){
	arrival_times = arrival_times_;
	Nt = param.Nt;
	t_min=param.t_min;
	t_max=param.t_max;
}

void Results::post_processing(){
	// time scale of arrival times
	double min_time = t_min;//get_min_value(arrival_times);
	//double min_time=get_min_value_pos(arrival_times);
	// January 10 2018 normalizing with the max arrival time rather then the max time of the system. We need to verify that it still produces correct results.
	//double max_time = MAX_TIME;//get_max_value(arrival_times);  //
	double max_time = t_max;//get_max_value(arrival_times);  //
	//double max_time=get_max_value(arrival_times);  
	/*print(arrival_times);
	cout << "min arrival time = " << min_time << ", max arrival time = " << max_time << endl;*/
	vector<double> time = define_scale(min_time,max_time,Nt,SCALE_LOG);
	// cumulative distribution computation  with change Jan 2016
	double current_time,current_time_next; int nb_part,nb_part_pdf;
	for (int i=0;i<Nt;i++){
		// shifting Originally it was: current_time = time[i]; to:
		if (time[i]>0.0){		
			current_time = time[i];
			nb_part = 0;
			for (map<int,double>::iterator it=arrival_times.begin();it!=arrival_times.end();it++){
				if (it->second<=current_time) nb_part++;
			}
			cum_dist_times[current_time] = nb_part;
		}	
	}
	// pdf needs a shift the CDF is adding one quicker then the PDF
	for (int i=0;i<Nt-1;i++){
		if (time[i]>0.0){
			current_time = time[i];
			current_time_next = time[i+1];
			nb_part_pdf = 0;
			for (map<int,double>::iterator it=arrival_times.begin();it!=arrival_times.end();it++){
				if (current_time<=it->second&&it->second<current_time_next) {nb_part_pdf++;}
				if (it->second==current_time_next&&current_time_next==max_time){nb_part_pdf++;}
			}
			//pdf_part[current_time] = nb_part_pdf;
			pdf_part[current_time] = (double)nb_part_pdf/(current_time_next-current_time);	// modified by Guofeng and Delphine (2022/07/20)
		}
	}
	pdf_part[time[Nt-1]] = 0;
}

void Results::writing(string code_path){
	int nb_part = arrival_times.size();
	string output_file=code_path+"/Output/cdf.txt";
	ofstream output1(output_file.c_str(), std::ofstream::out);
	if (nb_part==0){
		for (std::map<double,int>::iterator it=cum_dist_times.begin(); it!=cum_dist_times.end(); it++){
			output1 << it->first << "	0" << endl;
		}
		output1.close();

		output_file=code_path+"/Output/pdf.txt";
		ofstream output2_zero(output_file.c_str(), std::ofstream::out);
		for (std::map<double,double>::iterator it=pdf_part.begin(); it!=pdf_part.end(); it++){
			output2_zero << it->first << "	0" << endl;
		}
		output2_zero.close();
		return;
	}
	for (std::map<double,int>::iterator it=cum_dist_times.begin(); it!=cum_dist_times.end(); it++){
		output1 << it->first << "	"<< (double)it->second/nb_part << endl;
	}
	output1.close();

	output_file=code_path+"/Output/pdf.txt";
	ofstream output2(output_file.c_str(), std::ofstream::out);
	for (std::map<double,double>::iterator it=pdf_part.begin(); it!=pdf_part.end(); it++){
		output2 << it->first << "	"<< (double)it->second/nb_part << endl;
	}
	output2.close();

	double mean_data=0,var_data=0;
	for (std::map<int,double>::iterator it=arrival_times.begin(); it!=arrival_times.end(); it++){
		mean_data+=it->second;
	}
	mean_data=(double)mean_data/nb_part;

	for (std::map<int,double>::iterator it=arrival_times.begin(); it!=arrival_times.end(); it++){
		var_data+=(it->second-mean_data)*(it->second-mean_data);
	}
	var_data=(double)var_data/nb_part;
}

void Results::writing(string code_path,int i){
	int nb_part = arrival_times.size();
	stringstream ss; ss << i;
	string output_file=code_path+"/cdf/cdf"+ss.str()+".txt";
	ofstream output1(output_file.c_str(), std::ofstream::out);
	if (nb_part==0){
		for (std::map<double,int>::iterator it=cum_dist_times.begin(); it!=cum_dist_times.end(); it++){
			output1 << it->first << "	0" << endl;
		}
		output1.close();

		output_file=code_path+"/pdf/pdf"+ss.str()+".txt";
		ofstream output2_zero(output_file.c_str(), std::ofstream::out);
		for (std::map<double,double>::iterator it=pdf_part.begin(); it!=pdf_part.end(); it++){
			output2_zero << it->first << "	0" << endl;
		}
		output2_zero.close();
		return;
	}
	for (std::map<double,int>::iterator it=cum_dist_times.begin(); it!=cum_dist_times.end(); it++){
		output1 << it->first << "	"<< (double)it->second/nb_part << endl;
	}
	output1.close();

	output_file=code_path+"/pdf/pdf"+ss.str()+".txt";
	ofstream output2(output_file.c_str(), std::ofstream::out);
	for (std::map<double,double>::iterator it=pdf_part.begin(); it!=pdf_part.end(); it++){
		output2 << it->first << "	"<< (double)it->second/nb_part << endl;
	}
	output2.close();

	double mean_data=0,var_data=0;
	for (std::map<int,double>::iterator it=arrival_times.begin(); it!=arrival_times.end(); it++){
		mean_data+=it->second;
	}
	mean_data=(double)mean_data/nb_part;

	for (std::map<int,double>::iterator it=arrival_times.begin(); it!=arrival_times.end(); it++){
		var_data+=(it->second-mean_data)*(it->second-mean_data);
	}
	var_data=(double)var_data/nb_part;
}

pointcpp<double> ReadCoordinates(string str_line){
	string x_str,y_str;
	int pos1 = 0;
	pos1=str_line.find("(");
	int pos2=str_line.find(",",pos1);
	x_str=str_line.substr(pos1+1,pos2-1);
	pos1 = pos2;
	pos2 = str_line.find(")",pos1);
	y_str = str_line.substr(pos1+1,pos2-pos1-1);
	return pointcpp<double>(atof(x_str.c_str()),atof(y_str.c_str()));
}
