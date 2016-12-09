/*
 * Copyright © Alexandre Boulch, Renaud Marlet, Ecole Nationale des Ponts et Chaussees
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify,
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following
 * conditions:

 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.

 * The Software is provided "as is", without warranty of any kind, express or implied,
 * including but not limited to the warranties of merchantability, fitness for a
 * particular purpose and noninfringement. In no event shall the authors or copyright
 * holders be liable for any claim, damages or other liability, whether in an action of
 * contract, tort or otherwise, arising from, out of or in connection with the software
 * or the use or other dealings in the Software.
 */


#include <iostream>
#include <fstream>
#include "surface_comparator_MW.h"
#include "surface_comparator_KS.h"

#include <vector>
#include <time.h>


using namespace std;

class Point{
public:
	double x,y,z;

	Point():x(0), y(0), z(0){};
	Point(double _x, double _y, double _z):x(_x), y(_y), z(_z){}

	double dot(const Point& pt) const{
		return x*pt.x+y*pt.y+z*pt.z;
	}


	double norm(){
		return sqrt(x*x+y*y+z*z);
	}
};


Point operator+(const Point& pt1, const Point& pt2){
	return Point(pt1.x+pt2.x, pt1.y+pt2.y, pt1.z+pt2.z);
}

Point operator-(const Point& pt1, const Point& pt2){
	return Point(pt1.x-pt2.x, pt1.y-pt2.y, pt1.z-pt2.z);
}

Point operator*(double v, const Point& pt1){
	return Point(pt1.x*v, pt1.y*v, pt1.z*v);
}
Point operator*( const Point& pt1,double v){
	return v*pt1;
}

class Plane{
public:
	Point point;
	Point normal;
	
	double distance(const Point& p) const {
		return fabs(normal.dot(p-point));
	}
};

class Cylinder{
public:
	Point point;
	Point direction;
	double radius;
	
	double distance(const Point& p) const {
		Point v = p - point;
		return fabs(radius - (v - (v.dot(direction))*direction).norm());
	}
};

class Sphere{
public:
	Point point;
	double radius;
	
	double distance(const Point& p) const {
		return fabs(radius-(p-point).norm());
	}
};

typedef SurfaceComparatorMW<Point, Plane, Plane> SCMW_plane;
typedef SurfaceComparatorKS<Point, Plane, Plane> SCKS_plane;

typedef SCMW_plane::Point Point;
typedef SCMW_plane::PointSet PointSet;

void perturb_point_cloud(PointSet& P, double noise_position){
	for(int i=0; i<P.size(); i++){
		//box Muller method for normal variable
		double U = ((double)rand())/RAND_MAX ;
		double V = ((double)rand())/RAND_MAX ;
		double dir1 = sqrt(-2*log(U))*cos(2*3.14159265*V)*noise_position;
		U = ((double)rand())/RAND_MAX ;
		V = ((double)rand())/RAND_MAX ;
		double dir2 = sqrt(-2*log(U))*cos(2*3.14159265*V)*noise_position;
		U = ((double)rand())/RAND_MAX ;
		V = ((double)rand())/RAND_MAX ;
		double dir3 = sqrt(-2*log(U))*cos(2*3.14159265*V)*noise_position;
		P[i].x += dir1;
		P[i].y += dir2;
		P[i].z += dir3;
	}
}

void generate_horizontal_plane(PointSet & P,
		int nbr_points)
{
	assert(nbr_points>0);
	P.resize(nbr_points);
	for(int n=0; n<nbr_points; n++){
		double x,y;
		x = (rand()+0.f)/RAND_MAX*2.f - 1.f;
		y = (rand()+0.f)/RAND_MAX*2.f - 1.f;
		P[n] =Point(x,y,0);
	}
}


void translate_acquisition(PointSet& acq, const Point &translation_vector){
	for(int n=0; n<acq.size(); n++){
		acq[n] =acq[n]+ translation_vector;
	}
}


int main(){

	srand(time(NULL));
	cout << "Primitive fusion" << endl;
	
	PointSet P1,P2;

	double dist_plane =0.01;
	int nbr_pts = 1000;
	double noise = 0.03;
	
	for(int i=0; i<nbr_pts; i++){

		Point pt1((rand()+0.)/RAND_MAX, (rand()+0.)/RAND_MAX, 0);
		Point pt2((rand()+0.)/RAND_MAX, (rand()+0.)/RAND_MAX, dist_plane);
		
		//box Muller method for normal variable
		double U = ((double)rand())/RAND_MAX ;
		double V = ((double)rand())/RAND_MAX ;
		double dir1 = sqrt(-2*log(U))*cos(2*3.14159265*V)*noise;
		U = ((double)rand())/RAND_MAX ;
		V = ((double)rand())/RAND_MAX ;
		double dir2 = sqrt(-2*log(U))*cos(2*3.14159265*V)*noise;

		pt1.z += dir1;
		pt2.z += dir2;
		
		P1.push_back(pt1);
		P2.push_back(pt2);
		
	}
	
	
	Plane pl1,pl2;
	pl1.point = Point(0,0,0);
	pl2.point = Point(0,0,dist_plane);
	pl1.normal = pl2.normal = Point(0,0,1);
	

	cout << "SCMW" << endl;
	SCMW_plane scmw(&pl1,&P1,&pl2,&P2);
	cout << "MW: " <<scmw.apply_MW_test(0.01, 0.05) << endl;

	cout << "SCKS" << endl;
	SCKS_plane scks(&pl1,&P1,&pl2,&P2);
	cout << "KS: " << scks.apply_KS_test(0.01, 0.05) << endl;
}
