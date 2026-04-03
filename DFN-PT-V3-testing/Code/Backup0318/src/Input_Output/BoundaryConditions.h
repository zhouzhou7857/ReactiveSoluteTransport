/*
 * BoundaryConditions.h
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_

#include <string>
#include <map>
#include <vector>
#include "../Utilitaries/UblasStructures.h"

/************************************/
/* Functions for BoundaryConditions */
/************************************/
class BoundaryConditions{
public:
	std::string bc_type;
	double bc_value;
public:
	BoundaryConditions(){};
	BoundaryConditions(std::string bc_type_,double bc_value_=-1){bc_type=bc_type_;bc_value=bc_value_;}
	~BoundaryConditions(){};
	void print();
};

/************************************/
/* Functions for SourceTerms */
/************************************/
class SourceTerm{
public:
	double x_position;
	double y_position;
	double source_value;
public:
	SourceTerm(){};
	SourceTerm(double x_position_,double y_position_,double source_value_){x_position=x_position_;y_position=y_position_;source_value=source_value_;}
	~SourceTerm(){};
};

typedef std::map<int,BoundaryConditions> BoundaryConditionsDFN;
typedef std::map<std::string,BoundaryConditions> BoundaryConditionsDef;
typedef std::vector<SourceTerm> SourceTermsDef;
typedef std::map<double,BoundaryConditions> BoundaryConditionsMapCurv;


BoundaryConditionsDef ReadBoundaryConditions(std::string,SourceTermsDef&);
std::map<std::string,BoundaryConditionsDef> ReadBoundaryConditionsHydroPhys(std::string,SourceTermsDef&);

void print_BC(BoundaryConditionsDFN bc_map);
void print_BC(BoundaryConditionsMapCurv bc_map);

std::map<std::string,std::map<int,double> > ReturnBoundaryConditionsElec(BoundaryConditionsDef,ublas_matrix,ublas_matrix,std::map<int,double>,std::map<int,double>);



#endif /* BOUNDARYCONDITIONS_H_ */
