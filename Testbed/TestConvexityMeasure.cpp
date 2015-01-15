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
#include <ConvexityMeasure.h>
#include <Rank.h>

map<Triangulation::_Internal::_triangle*,int> 
getDepthMap(Triangulation::Triangulator& trmap)
{
	map<Triangulation::_Internal::_triangle*,int> dmap;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		Triangulation::_Internal::_triangle* tr = trmap.faces[i];
		if(tr->edges[0]->type==Triangulation::_Internal::Boundary || tr->edges[1]->type==Triangulation::_Internal::Boundary ||
			tr->edges[2]->type==Triangulation::_Internal::Boundary) 
		{
			dmap[trmap.faces[i]] = 0;
		}
		else
		{
			dmap[trmap.faces[i]] = trmap.edges.size(); //larger than possible depth
		}
	}
	while(true)
	{
		bool bChanged = false;
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int depth = dmap[trmap.faces[i]];
			for(int j=0; j<3; ++j)
			{
				Triangulation::_Internal::_edge* e = trmap.faces[i]->edges[j];
				Triangulation::_Internal::_triangle* f = e->faces[0]==trmap.faces[i] ? e->faces[1]: e->faces[0];
				if(f != NULL && dmap[f] + 1 < depth)
				{
					depth = dmap[f] + 1;
					bChanged = true;
				}
			}
			dmap[trmap.faces[i]] = depth;
		}
		if(bChanged == false) break;
	}
	return dmap;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [P E F] = TriangulationCarving(P)");
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
	float maxGap = -1.0;	
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(maxGap,prhs[1],classMode);
	} 

	Triangulation::Triangulator trmap(points);
	if(maxGap <= 0)
	{
		vector<float> vlength(trmap.edges.size());
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			vlength[i] = trmap.edges[i]->Length();
		}
		maxGap = Rank(vlength, 0.99);
		printf("99 percentile is %f\n", maxGap);
	}

	map<Triangulation::_Internal::_triangle*,int> depthMap = getDepthMap(trmap);
	int maxDepth = 0;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		int d = depthMap[trmap.faces[i]];
		if(d > maxDepth) maxDepth = d;
	}
	//printf("\nmax depth = %d, max int=%d\n", maxDepth, std::numeric_limits<int>::max());

	while(true)
	{
		bool bChanged = false;
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			Triangulation::_Internal::_edge* e = trmap.edges[i];
			if(e->type == Triangulation::_Internal::Boundary) // && fitness[i] < thres) 
			{
				if(e->Length() > maxGap)
				{
					trmap.RemoveAnyEdge(e);
					RemoveDanglingEdges(trmap);
					bChanged = true;
				}
			}
		}
		if(bChanged == false) break;
	}
	int maxDepth2 = 0;
	for(int i=0; i<trmap.deleted_faces.size(); ++i)
	{
		int d = depthMap[trmap.deleted_faces[i]];
		if(d > maxDepth2) maxDepth2 = d;
		//printf("%d: %d, %d\n", i, d, maxDepth2);
	}

	printf("%d vs. %d. Convexity measure = %f\n", maxDepth, maxDepth2, 1.0f - (float)maxDepth2/maxDepth);

	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trmap.faces.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[1]->p));
			int j3 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[2]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.edges.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], (int)trmap.edges[i]->type+1);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	mexUnlock();
}

