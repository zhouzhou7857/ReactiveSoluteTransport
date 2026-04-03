/*
 * Point_Cgal.cpp
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */


#include "Point_Cgal.h"
#include "Constantes.h"

using namespace std;


/**determine if two points are identic*/
int identic(const CgalPoint2D p1, const CgalPoint2D p2, const double eps){
	double x1 = p1.x(), x2 = p2.x(), y1 = p1.y(), y2 = p2.y();
	if((fabs(x1-x2)<eps)&&(fabs(y1-y2)<eps))  return 1; else return 0;
}

/**distance between two points*/
double distance_2D(const CgalPoint2D p1, const CgalPoint2D p2){
	double value = (p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y());
	return std::sqrt(value);
}

void print(CgalPoint2D pt){
	cout << "(" << pt.x() << "," << pt.y() << ")" << endl;
}

bool define(CgalPoint2D pt){
	if (pt.x()==NOT_DEFINED||pt.y()==NOT_DEFINED){return false;}
	return true;
}
