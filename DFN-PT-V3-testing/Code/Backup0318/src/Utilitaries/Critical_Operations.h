/*
 * Critical_Operations.h
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#ifndef CRITICAL_OPERATIONS_H_
#define CRITICAL_OPERATIONS_H_


/** compute the cross product of three points. */
inline double segment_cross_product( CgalPoint2D p1, CgalPoint2D p2, CgalPoint2D p3){
  return( ( p2.x() - p1.x()) * ( p3.y() - p1.y()) - ( p3.x() - p1.x()) * ( p2.y() - p1.y()));
}

/** test if a point belongs to a segment.*/
inline bool point_belongs_to_segment(const CgalSegment2D& seg, const CgalPoint2D& point, double EPS = 1e-10){
	CgalPoint2D pt1 = seg.source(), pt2 = seg.target();
	double xa = pt1.x();
	double ya = pt1.y();
	double xb = pt2.x();
	double yb = pt2.y();
	double xc = point.x();
	double yc = point.y();
	// vertical vector
	if(fabs(xc-xa)<EPS){
		// A==C
		if(fabs(yc-ya)<EPS) return true;
		else return (fabs(xb-xa)<EPS)&&((yc-ya)/(yb-ya)>-EPS)&&((yc-ya)/(yb-ya)<1.+EPS);
	}else if(fabs(yc-ya)<EPS){ // horizontal vector
		// A==C
		if(fabs(xc-xa)<EPS) return true; // useless (should be previously seen)
		else return (fabs(yb-ya)<EPS)&&((xc-xa)/(xb-xa)>-EPS)&&((xc-xa)/(xb-xa)<1.+EPS);
	}else{ // non vertical non horizontal vector
		// are the two vector colinear and C between A and B
		double k1 = (xc-xa)/(xb-xa);
		double k2 = (yc-ya)/(yb-ya);
		return (fabs(k1-k2)<EPS)&&(k1>-EPS)&&(k1<1.+EPS);
	}
}

inline bool Points_Equal( CgalPoint2D p1, CgalPoint2D p2, double diff = 1e-8){
	if((fabs(p1.x()-p2.x()) < diff) && (fabs(p1.y()-p2.y()) < diff))
		return true;
	else
		return false;
}



#endif /* CRITICAL_OPERATIONS_H_ */
