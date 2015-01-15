#ifndef ___LEVENBERG_MARQUARDT_H___
#define ___LEVENBERG_MARQUARDT_H___

#include <vector>
using namespace std;

#include <derivify.h>

bool
LevenbergMarquardt(const vector<double>& vvalues, //given values
				   vector<double>& vEstimates, //parameters to be derived
				   double (*func)(const vector<double>&, int n,  const void*), //function for evaluation
				   surreal (*funcD)(const vector<surreal>&, int n, const void*), //function for derivatives
				   int maxIter,
				   double weight,
				   const void* pParam = 0, //general purpose parameters for modeling
				   const double* lBounds = 0, //lower bounds of each estimate
				   const double* uBounds = 0,  //upper bounds of each estimate
				   const double* pSigma = 0,	//weights for each square term
				   ostream* out = 0
				   );

#endif /* ___LEVENBERG_MARQUARDT_H___ */