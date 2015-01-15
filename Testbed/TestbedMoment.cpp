#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <Triangulation.h>
#include <FragmentInfo.h>
#include <TriangulationHelper.h>
#include <DistanceMapUtility.h>
#include <TriangulationSaliencyUtility.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: Testbed(vals)");
		return;
	}

	vector<float> vals;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(vals, prhs[0], classIdP, ndimP, &dims);
	}

	vector<CParticleF> pnts;
	for(int i=0; i<dims[0]; i++)
	{
		float x = GetData2(vals, i, 0, dims[0], dims[1], 0.0f);
		float y = GetData2(vals, i, 1, dims[0], dims[1], 0.0f);
		float v = GetData2(vals, i, 2, dims[0], dims[1], 0.0f);
		pnts.push_back(CParticleF(x, y, 0, v));
	}
	vector<float> eg = Moment(pnts);
	for(int i=0; i<eg.size(); ++i)
	{
		printf("%f\n", eg[i]);
	}
	mexUnlock();
}

