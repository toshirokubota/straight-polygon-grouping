#ifndef ___NON_MAXIMUM_SUPPRESSION_H___
#define ___NON_MAXIMUM_SUPPRESSION_H___
#include <vector>
using namespace std;
#include <szMexUtility.h>
#include <szmexutilitytemplate.h>

vector<float>
CannyNonMaximumSuppression(const vector<float>& im,
                           const vector<float>& dy, 
						   const vector<float>& dx,
						   const int* dims);

vector<unsigned char>
CannyThinning(const vector<unsigned char>& image,
			  const int* dims);

vector<unsigned char>
thinning(const vector<float>& im, 
		 const vector<float>& dx,
		 const vector<float>& dy,
		 const int* dims);

#endif /* ___NON_MAXIMUM_SUPPRESSION_H___ */