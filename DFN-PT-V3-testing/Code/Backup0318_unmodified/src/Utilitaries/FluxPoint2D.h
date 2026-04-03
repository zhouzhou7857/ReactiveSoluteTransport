/*
 * FluxPoint2D.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef FLUXPOINT2D_H_
#define FLUXPOINT2D_H_

#include "Point_Cgal.h"


// mesh extremities definition
class FluxPoint2D{
public:
	CgalPoint2D p;
	int index;
public:
	FluxPoint2D(){};
	virtual ~FluxPoint2D(){};
	FluxPoint2D(CgalPoint2D,int index_=-1);
	FluxPoint2D(double,double,int index_=-1);
};


#endif /* FLUXPOINT2D_H_ */
