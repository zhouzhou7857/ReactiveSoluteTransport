/*
 * Point_Cgal.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef POINT_CGAL_H_
#define POINT_CGAL_H_

// CGAL arithmetic
#include <CGAL/Cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/intersection_2.h>

// CGAL Objects
#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Origin.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Vector_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Sphere_3.h>

//standards
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>


/** @addtogroup UTILITARY_GEOMETRY_CGAL UTILITARY GEOMETRY CGAL
 *  @{
 */

#define CGAL_DOUBLE

#ifdef CORE
typedef CORE::Expr real_NT;
#endif
#ifdef CGAL_DOUBLE
typedef double real_NT;
#endif

inline double convert_to_double( real_NT x){
#ifdef CORE
	return x.doubleValue();
#endif
#ifdef CGAL_DOUBLE
	return x;
#endif
}

typedef CGAL::Cartesian<real_NT>   K_CGAL;

#define CGAL_DISPLAY 0


/** A 2D point in CGAL.*/
typedef CGAL::Point_2<K_CGAL>                CgalPoint2D;
/** A 3D point in CGAL.*/
typedef CGAL::Point_3<K_CGAL>                CgalPoint3D;
/** A 2D triangle in CGAL.*/
typedef CGAL::Triangle_2<K_CGAL>                CgalTriangle2D;
/** A 2D polygon in CGAL.*/
typedef CGAL::Polygon_2<K_CGAL>				CgalPolygon2D;
/** A 2D polygon with holes in CGAL.*/
typedef CGAL::Polygon_with_holes_2<K_CGAL>			CgalPolygonWithHoles2D;
/** A 2D segment in CGAL.*/
typedef CGAL::Segment_2<K_CGAL>              CgalSegment2D;
/** A 3D plane in CGAL.*/
typedef CGAL::Plane_3<K_CGAL>				CgalPlane;
/** A 2D line in CGAL.*/
typedef CGAL::Line_2<K_CGAL>					CgalLine2D;
/** A 3D line in CGAL.*/
typedef CGAL::Line_3<K_CGAL>					CgalLine3D;
/** A 2D vector in CGAL.*/
typedef CGAL::Vector_2<K_CGAL>               CgalVector2D;
/** A 3D vector in CGAL.*/
typedef CGAL::Vector_3<K_CGAL>               CgalVector3D;
/** A 3D polyhedron in CGAL.*/
typedef CGAL::Polyhedron_3<K_CGAL>           CgalPolyhedron;
/** A 3D sphere in CGAL.*/
typedef CGAL::Sphere_3<K_CGAL>				CgalSphere;
/** A 3D affine transformation in CGAL.*/
typedef CGAL::Aff_transformation_3<K_CGAL>   Transformation3D;

/** computes the length of a Point2D for compatibility. */
inline real_NT length_(CgalPoint2D & s){
	real_NT x = 0;
	return x;
}

int identic(const CgalPoint2D,const CgalPoint2D,const double);
double distance_2D(const CgalPoint2D,const CgalPoint2D);
void print(CgalPoint2D);
bool define(CgalPoint2D);

#endif /* POINT_CGAL_H_ */
