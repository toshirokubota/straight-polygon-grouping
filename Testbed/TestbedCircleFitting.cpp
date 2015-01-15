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
#include <ContourEQW.h>
#include <szParticleF.h>
#include <TriangulationHelper.h>
#include <CircleFitting.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [] = Testbed(P)");
		return;
	}

	//Points
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}
	double radius=1.0, centerx=0.0, centery=0.0;
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(radius,prhs[1],classMode);
	}
	if(nrhs>=3) 
	{
		mxClassID classMode;
		ReadScalar(centerx,prhs[2],classMode);
	}
	if(nrhs>=4) 
	{
		mxClassID classMode;
		ReadScalar(centery,prhs[3],classMode);
	}
	int numiter=10;
	double rate=0.1;
	if(nrhs>=5) 
	{
		mxClassID classMode;
		ReadScalar(numiter,prhs[4],classMode);
	}
	if(nrhs>=6) 
	{
		mxClassID classMode;
		ReadScalar(rate,prhs[5],classMode);
	}

	printf("Initial value: Radius = %f, center = (%f, %f)\n", radius, centerx, centery);
	bool bSuccess = fitCircle(points, radius, centerx, centery, numiter, rate);
	double error = fittingCircleError(points, radius, centerx, centery);
	printf("Success = %d, Radius = %f, center = (%f, %f), error = %f\n", bSuccess, radius, centerx, centery, error);


	/*if(nlhs >= 1)
	{
		const int dims[] = {lmax.size(), 3};
		vector<float> F(dims[0]*dims[1], 0);
		for(int i=0; i<lmax.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], lmax[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], lmax[i].m_Y);
			SetData2(F, i, 2, dims[0], dims[1], lmax[i].m_Life);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StoreData(D, mxSINGLE_CLASS, 2, dimsD);
	}
	if(nlhs >= 3)
	{
		vector<float> C(dmap.size());
		CentricityTransform(C, dmap, dims);
		plhs[2] = StoreData(C, mxSINGLE_CLASS, 2, dims);
	}*/

	mexUnlock();
}

