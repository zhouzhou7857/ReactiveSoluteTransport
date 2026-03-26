/*
 * DFNComputation.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef DFN_COMPUTATION_H_
#define DFN_COMPUTATION_H_

#include "Domain.h"
#include "NetworkMeshes.h"
#include "../Input_Output/BoundaryConditions.h"

BoundaryConditionsDFN ReturnBoundCondDFN(Domain,BoundaryConditionsDef,BordersMap,NodesMap);
ublas_vector DFNComputationClassic(NetworkMeshes,BoundaryConditionsDFN);
ublas_vector DFNComputationClassic(NetworkMeshes,BoundaryConditionsDFN,SourceTermDFN,std::string option="conductivity",std::set<int> connected_nodes=std::set<int>());
ublas_vector ComputeFlowVelocities(NetworkMeshes&,Parameters,Domain,std::string option="neumann");


#endif /* DFN_COMPUTATION_H_ */
