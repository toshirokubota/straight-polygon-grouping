#ifndef ___CHAMFER_DISTANCE_H___
#define ___CHAMFER_DISTANCE_H___
#include <mex.h>
#include <vector>
using namespace std;

void
ChamferDistance(vector<float>& dmap,
				const vector<unsigned char>& fg,
				const int* dims);

void
CenterMaximalBall(vector<unsigned char>& center,
				  const vector<float>& dchamfer,
				  const vector<unsigned char>& fg,
				  const int* dims);


void
doChamferDistance(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[]);

#endif /* ___CHAMFER_DISTANCE_H___ */