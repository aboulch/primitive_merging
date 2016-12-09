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

#ifndef SURFACE_COMPARISON_KS_HEADER
#define SURFACE_COMPARISON_KS_HEADER

#include <vector>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <map>


template<class PtClass>
class MinimalClassDistKS{
public:
	double distance(const PtClass& p)const {return -1;}
};

template<class PtClass, class Surface1, class Surface2, class Surface3=MinimalClassDistKS<PtClass> >
class SurfaceComparatorKS{
	
	
public:	
	typedef  unsigned int uint;
	typedef PtClass Point;
	typedef std::vector<Point> PointSet;
	
	//vector for distributions
	std::vector<double> X,Y;
	
	
	//Surfaces and pointSets
	const Surface1* S1;
	const Surface2* S2;
	const Surface3* S3;
	const PointSet* pP1;
	const PointSet* pP2;
	
	SurfaceComparatorKS(const Surface1* s1,const PointSet* pp1,
						const Surface2* s2,const PointSet* pp2,
					 const Surface3* s3):S1(s1), S2(s2), S3(s3), pP1(pp1), pP2(pp2)
					 {}
	SurfaceComparatorKS(Surface1* s1,const PointSet* pp1, 
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
		
		if(typeid(S3).name() == typeid(const MinimalClassDistKS<Point>*).name()){//use surface 2
			generateXY_for_one_pointSet<Surface1,Surface2>(sample_size,*pP1,S1,S2);
			generateXY_for_one_pointSet<Surface2,Surface1>(sample_size,*pP2,S2,S1);
		}else{//use surface 3
			generateXY_for_one_pointSet<Surface1,Surface3>(sample_size,*pP1,S1,S3);
			generateXY_for_one_pointSet<Surface2,Surface3>(sample_size,*pP2,S2,S3);
		}
		
		std::sort(X.begin(), X.end());
		std::sort(Y.begin(), Y.end());
	}
		
	
	
	double Kolmogorov_Smirnov(){
		
		double Dmax = 0;
		double current_D=0;
		int current_rankX = 0;
		int current_rankY = 0;
		
		while(current_rankX < X.size() && current_rankY < Y.size()){
			if(X[current_rankX] < Y[current_rankY]){
				current_D += 1./X.size();
				current_rankX++;
			}else if(X[current_rankX] > Y[current_rankY]){
				current_D += -1./Y.size();
				current_rankY++;
			}else{
				//the two are equal
				if(rand()%2 == 0){//select direct
					current_D += 1./X.size();
					current_rankX++;
				}else{//select cross
					current_D += -1./Y.size();
					current_rankY++;
				}
			}
			if(fabs(current_D)>Dmax){
				Dmax = fabs(current_D);
			}
		}
		while(current_rankX < X.size()){
			current_D += 1./X.size();
			current_rankX++;
		}
		if(fabs(current_D)>Dmax){
			Dmax = fabs(current_D);
		}
		while(current_rankY < Y.size()){
			current_D += -1./Y.size();
			current_rankY++;
		}
		if(fabs(current_D)>Dmax){
			Dmax = fabs(current_D);
		}
		return Dmax/sqrt(2./X.size());
	}
	
	double KS_threshold(double alpha){
		double c_alpha=0;
		std::map<double,double> alpha_to_c_alpha;
		alpha_to_c_alpha[0.10] = 1.22;
		alpha_to_c_alpha[0.05] = 1.36;
		alpha_to_c_alpha[0.01] = 1.63;
		alpha_to_c_alpha[0.025] = 1.48;
		alpha_to_c_alpha[0.005] = 1.73;
		alpha_to_c_alpha[0.001] = 1.95;
		if(alpha_to_c_alpha.count(alpha)>0){
			c_alpha = alpha_to_c_alpha[alpha];
		}
		return c_alpha;
	}

	//apply test
	bool apply_KS_test(
		const double alpha, 
		const double epsilon){
		
		//rand init
		srand(time(NULL));
		
		//get the size of the samples
		uint nbr_points = sample_size(alpha,epsilon);
		
		//create vector X and Y
		generateXY(nbr_points);
		
		double c_alpha=KS_threshold(alpha);
		return Kolmogorov_Smirnov() < c_alpha;
	}
	
	//get statistic
	double get_KS(
		const double alpha, 
		const double epsilon){
		
		//rand init
		srand(time(NULL));
		
		//get the size of the samples
		uint nbr_points = sample_size(alpha,epsilon);
		
		//create vector X and Y
		generateXY(nbr_points);
		
		return Kolmogorov_Smirnov();
	}
};

#endif
