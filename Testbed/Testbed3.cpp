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

void
CentricityTransform(vector<float>& C, 
					const vector<float>& D,
					const int* dims)
{
	int xoff[] = {-1, 0, 1, -1, 1, -1, 0, 1};
	int yoff[] = {-1, -1, -1, 0, 0, 1, 1, 1};
	for(int i=0; i<dims[1]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			int larger = 0;
			int same = 0;
			float dval = GetData2(D, j, i, dims[0], dims[1], (float)0);
			for(int k=0; k<8; ++k)
			{
				float dval2 = GetData2(D, j+xoff[k], i+yoff[k], dims[0], dims[1], (float)0);
				if(dval > dval2) larger++;
				else if(dval==dval2) same++;
			}
			if(same <= 2) 
			{
				larger += same; //permissible
			}
			SetData2(C, j, i, dims[0], dims[1], (float)larger/8);
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 2 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [L D] = Testbed(dmap, P)");
		return;
	}

	//Points
	vector<float> dmap;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(dmap, prhs[0], classIdP, ndimP, &dims);
	}
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}

	vector<CParticleF> lmax = localMaximaPoints(dmap, 3.0f, dims); 
	const int dimsD[] = {points.size(), lmax.size()};
	vector<float> D(dimsD[0]*dimsD[1], 0);
	for(int i=0; i<points.size(); ++i)
	{
		for(int j=0; j<lmax.size(); ++j)
		{
			float d = Distance(points[i], lmax[j]) - lmax[j].m_Life;
			SetData2(D, i, j, dimsD[0], dimsD[1], d);
		}
	}

	if(nlhs >= 1)
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
	}

	mexUnlock();
}

