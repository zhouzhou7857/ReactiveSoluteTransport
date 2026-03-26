/*
 * Pointcpp.h
 *
 *  Created on: Jun 27, 2012
 *      Author: delphine
 */

#ifndef POINTCPP_H_
#define POINTCPP_H_


#include "Point_Cgal.h"


template <class T>
class pointcpp {
public:
	T i;
	T j;
	T k;
public:
	pointcpp(T a=0, T b=0, T c=0):i(a),j(b),k(c){}
	~pointcpp(){};
	pointcpp(CgalPoint2D pt){i=pt.x();j=pt.y();k=0;}
	void define(T a=0, T b=0, T c=0){ i = a;j = b;k = c;}
	CgalPoint2D CgalPoint(){return CgalPoint2D(i,j);}
	void print(){std::cout << "(" << i << "," << j  << ")" << std::endl;}
	pointcpp operator+(const pointcpp<T> & P)const{return pointcpp<T>(i + P.i, j + P.j, k + P.k);}
	pointcpp operator*(T temp)const{return pointcpp<T>(i*temp, j*temp, k*temp); }
	bool operator<(pointcpp<T> & P)const{
		if(i < P.i)
			return true;
		else if(i==P.i)
			if(j < P.j)
				return true;
			else if(j==P.j)
				if(k < P.k)
					return true;
		return false;
	}

	bool operator>(const pointcpp<T> & P)const{
		return P<(*this);
	}

	T Point_Distance(pointcpp<T> P)const{
		T temp = (i-P.i)*(i-P.i)+(j-P.j)*(j-P.j)+(k-P.k)*(k-P.k);
		temp = sqrt(temp);
		return temp;
	}
	bool identic(pointcpp<T> p, double eps)const{
		return((fabs(p.i-i)<eps)&&(fabs(p.j-j)<eps)&&(fabs(p.k-k)<eps));
	}
};

template< typename T>
pointcpp<T> operator*(const T temp, const pointcpp<T> & P){return P*temp;}

template <typename T>
inline bool operator < ( const pointcpp<T> & p1, const pointcpp<T> & p2) {
	if(p1.i!=p2.i) {
		if(p1.i<p2.i) return 1;
		else		return 0;}
	if(p1.j!=p2.j) {
		if(p1.j<p2.j) return 1;
		else		return 0;}
	if(p1.k!=p2.k) {
		if(p1.k<p2.k) return 1;
		else		return 0;}
	return 0;
}

template <typename T>
inline bool operator < ( pointcpp<T> & p1, pointcpp<T> & p2) {
	if(p1.i!=p2.i) {
		if(p1.i<p2.i) return 1;
		else		return 0;}
	if(p1.j!=p2.j) {
		if(p1.j<p2.j) return 1;
		else		return 0;}
	if(p1.k!=p2.k) {
		if(p1.k<p2.k) return 1;
		else		return 0;}
	return 0;
}



#endif /* POINTCPP_H_ */
