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

#ifndef SURFACE_COMPARATOR_MW_HEADER
#define SURFACE_COMPARATOR_MW_HEADER

#include <vector>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <typeinfo>

#include <boost/math/distributions/normal.hpp>

template<class PtClass>
class MinimalClassDistMW{
public:
	double distance(const PtClass& p)const {return -1;}
};

template<class PtClass, class Surface1, class Surface2, class Surface3=MinimalClassDistMW<PtClass> >
class SurfaceComparatorMW{
	
public:
	
	typedef  unsigned int uint;
	typedef PtClass Point;
	typedef typename std::vector<Point> PointSet;
	
	//vector for distributions
	std::vector<double> X,Y;
	
	//Surfaces and pointSets
	const Surface1* S1;
	const Surface2* S2;
	const Surface3* S3;
	const PointSet* pP1;
	const PointSet* pP2;
	
	SurfaceComparatorMW(const Surface1* s1,const PointSet* pp1,
						const Surface2* s2,const PointSet* pp2,
					 const Surface3* s3):S1(s1), S2(s2), S3(s3), pP1(pp1), pP2(pp2)
					 {}
	SurfaceComparatorMW(Surface1* s1,const PointSet* pp1, 
						Surface2* s2,const PointSet* pp2):S1(s1), S2(s2), pP1(pp1), pP2(pp2){}
	
	
	uint sample_size(double alpha, double epsilon){
		//from Dvoretsy-Keifer-Wolfowitz
		return ceil(log(2./alpha))/(2*epsilon*epsilon);
	}
	
	inline uint minUint(const uint& i1, const uint& i2){
		if(i1 < i2)	return i1;
		return i2;
	}
	
	template<class SurfaceP, class SurfaceOther>
	void generateXY_for_one_pointSet(uint sample_size,
		const PointSet& P,
		const SurfaceP* SP,
		const SurfaceOther* SOther){
		//have to divide sample size by 2 for the two set of points
		uint nbr_points_prim = ceil(sample_size/2.);
		
		if(P.size() <= nbr_points_prim){
			//less points than desired => use all points
			for(uint i=0; i<P.size(); i++){
				X.push_back(SP->distance(P[i]));
				Y.push_back(SOther->distance(P[i]));
			}
		}else{
			//shuffle a copy of P1
			//slow and memory greedy, should be fixed
			PointSet copyP = P;
			for(uint i=0; i<nbr_points_prim; i++){
				uint idSwap =i+rand()%(P.size()-i);
				Point temp = copyP[idSwap];
				copyP[idSwap] = copyP[i];
				copyP[i] = temp;
			}
			for(uint i=0; i<nbr_points_prim; i++){
				X.push_back(SP->distance(copyP[i]));
				Y.push_back(SOther->distance(copyP[i]));
			}
		}
	}
	
	void generateXY(uint sample_size){
		//clear vector
		X.clear();
		Y.clear();
		
		if(typeid(S3).name() == typeid(const MinimalClassDistMW<Point>*).name()){//use surface 2
			generateXY_for_one_pointSet<Surface1,Surface2>(sample_size,*pP1,S1,S2);
			generateXY_for_one_pointSet<Surface2,Surface1>(sample_size,*pP2,S2,S1);
		}else{//use surface 3
			generateXY_for_one_pointSet<Surface1,Surface3>(sample_size,*pP1,S1,S3);
			generateXY_for_one_pointSet<Surface2,Surface3>(sample_size,*pP2,S2,S3);
		}
		
		std::sort(X.begin(), X.end());
		std::sort(Y.begin(), Y.end());
	}
	
	double Mann_Whitney(){

		//compute x rank
		double total_rankX = 0;
		double total_rankY = 0;
		int current_rankX = 0;
		int current_rankY = 0;
		while(current_rankX < X.size() && current_rankY < Y.size()){
			if(X[current_rankX] < Y[current_rankY]){
				current_rankX++;
				total_rankX+=current_rankX+current_rankY;
				continue;
			}
			if(X[current_rankX] > Y[current_rankY]){
				current_rankY++;
				total_rankY+=current_rankX+current_rankY;
				continue;
			}
			total_rankX+=current_rankX+current_rankY+0.5;
			total_rankY+=current_rankX+current_rankY+0.5;
			current_rankX++;
			current_rankY++;
		}
		while(current_rankX < X.size()){
			current_rankX++;
			total_rankX+=current_rankX+current_rankY;
		}
		while(current_rankY < X.size()){
			current_rankY++;
			total_rankY+=current_rankX+current_rankY;
		}

		double Uy = double(total_rankY - (Y.size())*(Y.size()+1)/2);
		double Xi_X_Y = (Uy-(X.size()*Y.size()/2.))/sqrt(X.size()*Y.size()*(X.size()+Y.size()+1)/12.);
		

		return 1-boost::math::cdf(boost::math::normal(0,1), fabs(Xi_X_Y));

	}
	
	bool apply_MW_test(
		const double alpha, 
		const double epsilon){
		
		//rand init
		srand(time(NULL));
		
		//get the size of the samples
		uint nbr_points = sample_size(alpha,epsilon);
		
		//create vector X and Y
		generateXY(nbr_points);
		
		return Mann_Whitney() > alpha;
	}
	
	double get_MW(
		const double alpha, 
		const double epsilon){
		
		//rand init
		srand(time(NULL));
		
		//get the size of the samples
		uint nbr_points = sample_size(alpha,epsilon);
		//create vector X and Y
		generateXY(nbr_points);
		
		
		return Mann_Whitney();
	}
};


#endif
