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
#include <map>
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
#include <DotGroupUtility.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Kruskal.h>
#include <cassert>

mxArray*
StoreContours(const vector<vector<CParticleF>>& polygons)
{
	const int dims[] = { polygons.size() };
	mxArray* cell = mxCreateCellArray(1, (mwSize*)dims);
	for (int i = 0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size();
		const int dimsC[] = { n, 3 };
		mxArray* ar = mxCreateNumericArray(2, (mwSize*)dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*)mxGetData(ar);
		for (int j = 0; j<n; ++j)
		{
			p[j] = polygons[i][j].m_X;
			p[n + j] = polygons[i][j].m_Y;
			p[2 * n + j] = polygons[i][j].m_Z;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

void
splitNow(vector<CParticleF>& points, int beg, int end,
		vector<int>& match,
		vector<vector<CParticleF>>& polygons)
{
	if (beg >= end) return;
	vector<CParticleF> poly;
	for (int i = beg; i <= end; )
	{
		poly.push_back(points[i]);
		if (match[i] > i)
		{
			splitNow(points, i+1, match[i], match, polygons);
			i = match[i]+1;
		}
		else
		{
			i++;
		}
	}
	polygons.push_back(poly);
}

vector<vector<CParticleF>> 
split(vector<CParticleF>& points)
{
	vector<vector<CParticleF>> polygons;
	vector<int> match(points.size(), -1);
	float eps = 1.0e-5;
	for (int i = 0; i < points.size(); ++i)
	{
		int d = points.size();
		for (int j = 0; j < points.size(); ++j)
		{
			if (i == j) continue;
			if (Distance(points[i], points[j]) < eps)
			{
				int k = Abs(i - j);
				if (k < d)
				{
					match[i] = j;
				}
			}
		}
	}

	splitNow(points, 0, points.size() - 1, match, polygons);

	return polygons;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "ClosedPolygons", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [P G] = ClosedPolygons(P)");
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

	vector<vector<CParticleF>> polygons = split(points);
	vector<vector<CParticleF>> keep;
	float eps = 1.0e-5;
	for (int i = 0; i < polygons.size(); ++i)
	{
		if (polygonArea(polygons[i]) < eps) continue;

		keep.push_back(polygons[i]);
	}

	if(nlhs >= 1)
	{
		const int dims[] = {points.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreContours(keep);
	}
	mexUnlock();
}

