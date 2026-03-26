/*
 * Segment.h
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#ifndef SEGMENT_H_
#define SEGMENT_H_

#include "Point_Cgal.h"



class Segment2D : public CgalSegment2D{
public:
	double length; ///< the length of the segment (used for inter_system)
	CgalPoint2D center; ///< the center of the segment (used for inter_system)
public:
	Segment2D();
	Segment2D(CgalPoint2D,double,double);
	Segment2D(CgalPoint2D p1, CgalPoint2D p2);
	Segment2D& operator=(const Segment2D& other);
	Segment2D(const Segment2D& other);
	double get_orientation() const;
	bool advanced_reject_intersection(const Segment2D&);
	bool immediate_reject_intersection(const Segment2D&);
	bool intersection_possible(const Segment2D&);
	bool intersection(const Segment2D&,CgalPoint2D&);
	double distance(const CgalPoint2D&);
	bool IdenticExtremities(CgalPoint2D);
};

bool Intersect(Segment2D,Segment2D,CgalPoint2D&);
bool Intersect(CgalLine2D,Segment2D,CgalPoint2D&);
bool Intersect(CgalLine2D,CgalLine2D,CgalPoint2D&);
void print(Segment2D);
void print(CgalLine2D);

#endif /* SEGMENT_H_ */
