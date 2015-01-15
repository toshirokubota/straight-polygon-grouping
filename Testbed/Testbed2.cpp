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

float
HausdorffDistance(const vector<CParticleF>& path1, const vector<CParticleF>& path2)
{
	float distance = 0;
	for(int i=0; i<path1.size(); ++i)
	{
		float mind = std::numeric_limits<float>::infinity();
		for(int j=0; j<path2.size(); ++j)
		{
			float d = Distance(path1[i], path2[j]);
			mind = Min(d, mind);
		}
		distance = Max(mind, distance);
	}
	for(int i=0; i<path2.size(); ++i)
	{
		float mind = std::numeric_limits<float>::infinity();
		for(int j=0; j<path1.size(); ++j)
		{
			float d = Distance(path2[i], path1[j]);
			mind = Min(d, mind);
		}
		distance = Max(mind, distance);
	}
	return distance;
}

float
getEdgeCost(vector<CParticleF> path1,
			vector<CParticleF> path2)
{
	//vector<CParticleF> path1 = traceUpwardRepresentative(e->vertices[0]->p, dmap, dims);
	//vector<CParticleF> path2 = traceUpwardRepresentative(e->vertices[1]->p, dmap, dims);

	float c1 = HausdorffDistance(path1, path2);
	//return (Distance(path1[0], path2[0]) + Distance(path1[path1.size()-1], path2[path2.size()-1]))/2.0f;
	float c2 = Max(Distance(path1[0], path2[0]), Distance(path1[path1.size()-1], path2[path2.size()-1]));
	return Max(c1, c2);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: T = TriangleClustering(P, T)");
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

	vector<vector<CParticleF>> traces(points.size());
	for(int i=0; i<points.size(); ++i)
	{
		points[i].m_Life = GetData2(dmap, points[i].m_X, points[i].m_Y, dims[0], dims[1], 0.f);
		printf("Tracing from (%f,%f) with %f.\n", points[i].m_X, points[i].m_Y, points[i].m_Life);
		traces[i] = traceUpwardVariableNeighbor(points[i], dmap, 4, dims);
	}

	vector<ContourEQW> contours(points.size());
	for(int i=0; i<points.size(); ++i)
	{
		ContourEQW contour(traces[i].size());
		for(int j=0; j<traces[i].size(); ++j)
		{
			contour.X[j] = traces[i][j].m_X; //one based index
			contour.Y[j] = traces[i][j].m_Y; //one based index
			contour.Strength[j] = traces[i][j].m_Life; 
		}
		contours[i] = contour;
	}

	vector<float> cost(points.size()*points.size(), 0);
	for(int i=0; i<points.size(); ++i)
	{
		for(int j=i; j<points.size(); ++j)
		{
			float val = getEdgeCost(traces[i], traces[j]);
			SetData2(cost, i, j, points.size(), points.size(), val);
			SetData2(cost, j, i, points.size(), points.size(), val);
		}
	}

	if(nlhs >= 1)
	{
		plhs[0] = StoreContoursEQW(contours);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {points.size(), points.size()};
		plhs[1] = StoreData(cost, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

