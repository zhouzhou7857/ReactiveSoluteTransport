/*
 * FluxPoint2D.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */



#include "FluxPoint2D.h"



FluxPoint2D::FluxPoint2D(CgalPoint2D p_,int index_){
	p = p_;
	index = index_;
}
FluxPoint2D::FluxPoint2D(double x_,double y_,int index_){
	p = CgalPoint2D(x_,y_);
	index = index_;
}

