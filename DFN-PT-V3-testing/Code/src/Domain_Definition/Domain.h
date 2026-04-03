/*
 * Domain.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef DOMAIN_H_
#define DOMAIN_H_

#include "../Utilitaries/Pointcpp.h"
#include "../Utilitaries/Segment.h"

class Domain {
public:
	pointcpp<double> min_pt;
	pointcpp<double> max_pt;
public:
	Domain();
	virtual ~Domain();
	Domain(pointcpp<double>,pointcpp<double>);
	Domain(double,double);
	bool on_input_limit(pointcpp<double>);
	bool on_output_limit(pointcpp<double>);
	pointcpp<double> return_middle_input_limit();
	double domain_size_x(){return max_pt.i-min_pt.i;}
	double domain_size_y(){return max_pt.j-min_pt.j;}
	void get_extremities(pointcpp<double> & pmin, pointcpp<double> & pmax){pmin=min_pt; pmax=max_pt;}
	std::string ReturnBorder(CgalPoint2D);
	Segment2D ReturnBorder(std::string);
	double ReturnBorderCoordinate(CgalPoint2D,double&,std::string&,std::string&);
	bool IsInDomain(CgalPoint2D);
	bool IsInDomain(Segment2D);
	bool IntersectionBorders(Segment2D,CgalPoint2D&,CgalPoint2D&);
	bool IntersectionBorders(CgalLine2D,CgalPoint2D&,CgalPoint2D&);
	bool IntersectionSingleBorderExtended(CgalLine2D,std::string,CgalPoint2D&);
	bool IntersectionBordersExtended(CgalLine2D,CgalPoint2D&,CgalPoint2D&);
	bool SegmentIntersectDomain(Segment2D,Segment2D&);
	CgalPoint2D ReturnSymmetric(CgalPoint2D,CgalPoint2D);
};

#endif /* DOMAIN_H_ */
