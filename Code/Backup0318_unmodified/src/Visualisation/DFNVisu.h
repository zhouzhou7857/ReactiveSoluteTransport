/*
 * DFNVisu.h
 *
 *  Created on: Jul 24, 2013
 *      Author: roubinet
 */

#ifndef DFNVISU_H_
#define DFNVISU_H_

#include "../Domain_Definition/NetworkMeshes.h"


class DFNVisu{
public:
	NetworkMeshes net_mesh;
public:
	DFNVisu();
	~DFNVisu(){};
	DFNVisu(NetworkMeshes);
	void display_DFN();
};

class DFNVisuSeg{
public:
	MapSegmentTimesCum seg_results;
	MapSegmentVelocities segment_velocities;
	Domain domain;
	double snapshot_time;
public:
	DFNVisuSeg();
	~DFNVisuSeg(){};
	DFNVisuSeg(MapSegmentTimesCum,MapSegmentVelocities,Domain,double);
	void display_DFN_Seg();
};

void SetColor(double);
void SetColorBis(double);
void SetColor(double,double,double);



#endif /* DFNVISU_H_ */
