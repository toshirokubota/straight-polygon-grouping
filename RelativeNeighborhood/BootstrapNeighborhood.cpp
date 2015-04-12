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
#include <set>
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
#include <Rank.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "RelativeNeighborhood", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [P E] = RelativeNeighborhood(P)");
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
	float scale = 2.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(scale, prhs[1], classMode);
	}

	Triangulation::Triangulator trmap(points);
	
	map<Triangulation::_Internal::_vertex*, int> pmap;
	for (int i = 0; i<trmap.points.size(); ++i)
	{
		pmap[trmap.points[i]] = i;
	}
	map<Triangulation::_Internal::_edge*, int> emap;
	for (int i = 0; i<trmap.edges.size(); ++i)
	{
		emap[trmap.edges[i]] = i;
	}

	//find the scale (length threshold) based on the length of triangulated graph
	vector<pair<float, Triangulation::_Internal::_edge*>> pedges(trmap.edges.size());
	for (int i = 0; i < trmap.edges.size(); ++i)
	{
		pedges[i].first = trmap.edges[i]->Length();
		pedges[i].second = trmap.edges[i];
	}
	sort(pedges.begin(), pedges.end());
	float thres = pedges[points.size()].first;
	vector<bool>  bkeep(pedges.size());
	vector<bool>  bcovered(points.size(), 0);
	for (int i = 0; i < pedges.size(); ++i) 
	{
		if (pedges[i].first > thres) break;
		Triangulation::_Internal::_edge* e = pedges[i].second;
		for (int j = 0; j < 2; ++j) 
		{
			int k = pmap[e->vertices[j]];
			if (bcovered[k] == false)
			{
				bkeep[emap[e]] = true;
				bcovered[k] = true;
			}
		}
	}
		
	if (nlhs >= 1)
	{
		const int dims[] = { trmap.points.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { trmap.edges.size(), 3 };
		vector<int> F(dims[0] * dims[1], 0);
		for (int i = 0; i<trmap.edges.size(); ++i)
		{
			int j1 = pmap[trmap.edges[i]->vertices[0]];
			int j2 = pmap[trmap.edges[i]->vertices[1]];
			SetData2(F, i, 0, dims[0], dims[1], (j1 + 1));
			SetData2(F, i, 1, dims[0], dims[1], (j2 + 1));
			SetData2(F, i, 2, dims[0], dims[1], bkeep[i] ? 1 : 0);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	mexUnlock();
}

