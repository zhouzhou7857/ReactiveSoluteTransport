/*
 * Segment.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */


#include "Segment.h"
#include "Constantes.h"
#include "Critical_Operations.h"

Segment2D::Segment2D() : CgalSegment2D(){
	length = 0;
	center = CGAL::ORIGIN;
}

Segment2D::Segment2D(CgalPoint2D center_,double length_,double angle){
	length=length_;
	center=center_;
	CgalPoint2D p1(center.x()-(length*cos(angle)/2),center.y()-(length*sin(angle)/2));
	CgalPoint2D p2(center.x()+(length*cos(angle)/2),center.y()+(length*sin(angle)/2));
	*this = Segment2D(p1,p2);
}

Segment2D::Segment2D(CgalPoint2D p1, CgalPoint2D p2) : CgalSegment2D(p1,p2){
	length = CGAL::sqrt(CGAL::squared_distance(p1,p2));
	center = CGAL::midpoint(p1,p2);
}

Segment2D& Segment2D::operator =(const Segment2D& other){
	CgalSegment2D::operator=(other);
	length = other.length;
	center = other.center;
	return *this;
}

Segment2D::Segment2D(const Segment2D& other){
	CgalSegment2D::operator=(other);
	length = other.length;
	center = other.center;
}

double Segment2D::get_orientation() const {
	// use Pythagore cos(phi) = AB/AC with AC hypothenuse
	CgalPoint2D pt1 = source(), pt2 = target();
	double dist = CGAL::sqrt(CGAL::squared_distance(pt1,pt2)); // AC
	if(dist==0.) return 0.; // case A==C
	double x1 = pt1.x(), x2 = pt2.x();
	double cosangle = CGAL::abs(x1-x2)/dist; // phi
	if(cosangle<-1||cosangle>1) return 0.; // case error (not a triangle, normally never happens)
	double angle = acos(cosangle);
	double y1 = pt1.y(), y2 = pt2.y();
	if(((x1<x2)&&(y1<y2))||((x1>x2)&&(y1>y2)))
		return angle;
	else
		return 2*PI-angle;
}

bool Segment2D::advanced_reject_intersection(const Segment2D & Seg){
	double eps = 1e-12;
	CgalPoint2D pt1=source(), pt2=target(), pt3=Seg.source(), pt4=Seg.target();
	if(((segment_cross_product(pt1,pt2,pt3) * segment_cross_product(pt1,pt2,pt4)) < eps)
      &&((segment_cross_product(pt3,pt4,pt1) * segment_cross_product(pt3,pt4,pt2)) < eps))
		return( false);
	else return( true);
}

bool Segment2D::immediate_reject_intersection(const Segment2D & Seg){
	CgalPoint2D p1 = this->source();
	CgalPoint2D p2 = this->target();
	CgalPoint2D p3 = Seg.source();
	CgalPoint2D p4 = Seg.target();
	double eps = 1e-12;
	 if( ((((p4.x()-p1.x())>eps)&&((p4.x()-p2.x())>eps)&&
          ((p3.x()-p1.x())>eps)&&((p3.x()-p2.x()))>eps) ||
          (((p1.x()-p4.x())>eps)&&((p2.x()-p4.x())>eps)&&
          ((p1.x()-p3.x())>eps)&&((p2.x()-p3.x())>eps)))
          ||
          ((((p4.y()-p1.y())>eps)&&((p4.y()-p2.y())>eps)&&
          ((p3.y()-p1.y())>eps)&&((p3.y()-p2.y())>eps)) ||
          (((p1.y()-p4.y())>eps)&&((p2.y()-p4.y())>eps)&&
          ((p1.y()-p3.y())>eps)&&((p2.y()-p3.y())>eps))) )
          return( true);
  else return( false);
}

bool Segment2D::intersection_possible(const Segment2D & Seg){
	if(!immediate_reject_intersection(Seg))
		if(!advanced_reject_intersection(Seg))
			return true;
	return false;
}

bool Segment2D::intersection(const Segment2D & Seg, CgalPoint2D & s_inter){
	// test if the segments are connected
	bool test = intersection_possible(Seg);
	if(test){
		CGAL::Object result = CGAL::intersection(*this,Seg);
		CgalPoint2D point;
		CgalSegment2D segment;
		if (assign(point, result)){
			s_inter = point;
		}else if (assign(segment,result)){
			#ifdef DEBUG
			if(*this==Seg){
				std::cout << "warning: two identical fractures in the network\n";
			}else{
				std::cout << "warning: two colinear fractures\n";
			}
			#endif
			s_inter = segment[0];
		}
		else{
			// the segment are not connected but might be parallel and very near
			if(CGAL::squared_distance(*this,Seg) <EPSILON*EPSILON){
				CgalPoint2D seg_source = Seg.source();
				CgalPoint2D seg_target = Seg.target();
				if (this->distance(seg_source)<this->distance(seg_target))
					s_inter = seg_source;
				else
					s_inter = seg_target;
				return true;
			}else{
				return false;
			}
		}
	}
	return test;
}

// Returns either:
// - the distance to the segment if the projected point is in the segment
// - the distance to the closest segment extremity if the projected point is outside of the segment
double Segment2D::distance(const CgalPoint2D & p){
	CgalLine2D l( CgalSegment2D(*this));
	CgalPoint2D p_proj = l.projection (p);
	double d_p_proj = CGAL::sqrt(CGAL::squared_distance(p, p_proj));
	if(!point_belongs_to_segment(*this,p_proj)) {
		CgalPoint2D pt1 = source(), pt2 = target();
		double a = CGAL::sqrt(CGAL::squared_distance(pt1, p_proj));
		double b = CGAL::sqrt(CGAL::squared_distance(pt2, p_proj));
		if(a>b) return sqrt((pow(b,2)+pow(d_p_proj, 2))); else return sqrt((pow(a,2)+pow(d_p_proj, 2)));
	} else
		return d_p_proj;
}

bool Segment2D::IdenticExtremities(CgalPoint2D pt){
	if (Points_Equal(this->source(),pt)||Points_Equal(this->target(),pt)){return true;}
	return false;
}

void print(Segment2D seg){
	std::cout << "(" << seg.source().x() << "," << seg.source().y() << ") -> (" << seg.target().x() << "," << seg.target().y() << ")" << std::endl;
}

void print(CgalLine2D line){
	std::cout << "a = " << line.a() << ", b = " << line.b() << ", c = " << line.c() << std::endl;
}

bool Intersect(Segment2D seg1,Segment2D seg2,CgalPoint2D& inter){
	CGAL::Object result=CGAL::intersection(seg1,seg2);
	if (CGAL::assign(inter,result)){return true;}
	return false;
}

bool Intersect(CgalLine2D line,Segment2D seg,CgalPoint2D& inter){
	CGAL::Object result=CGAL::intersection(line,seg);
	if (CGAL::assign(inter,result)){return true;}
	return false;
}

bool Intersect(CgalLine2D line1,CgalLine2D line2,CgalPoint2D& inter){
	CGAL::Object result=CGAL::intersection(line1,line2);
	if (CGAL::assign(inter,result)){return true;}
	return false;
}


