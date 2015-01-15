#ifndef ___SZ_GENERALIZED_DISTANCE_TRANSFORM_H___
#define ___SZ_GENERALIZED_DISTANCE_TRANSFORM_H___

#include <vector>
using namespace std;

void
GeneralizedDistanceTransform(vector<double>& D, 
	                        const vector<double>& L, 
							double maxValue, //maximum valued allowed in L
							int ndim, 
							const int* dims);

void
GeneralizedDistanceTransformF(vector<float>& D, 
	                        const vector<float>& L, 
							float maxValue, //maximum valued allowed in L
							int ndim, 
							const int* dims);

#endif /*  ___SZ_GENERALIZED_DISTANCE_TRANSFORM_H___ */
