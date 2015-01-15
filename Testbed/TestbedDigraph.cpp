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

Triangulation::_Internal::_vertex*
oppositeVertex(Triangulation::_Internal::_triangle* face,
			   Triangulation::_Internal::_edge* edge)
{
	for(int i=0; i<3; ++i)
	{
		Triangulation::_Internal::_vertex* u = face->vertices[i];
		if(u != edge->vertices[0] && u != edge->vertices[1])
		{
			return u;
		}
	}
	return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [] = Testbed(C)");
		return;
	}

	vector<ContourEQW> contours;
	LoadContourEQW(contours, prhs[0]);

	int numSamplePoints = 5;	
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(numSamplePoints,prhs[1],classMode);
	} 

	vector<vector<FragmentInfo>> vvinfo = TriangulationPointsFixedInterval(contours, numSamplePoints);
	vector<CParticleF> points;
	vector<FragmentInfo> vinfo;
	for(int i=0; i<vvinfo.size(); ++i)
	{
		for(int j=0; j<vvinfo[i].size(); ++j)
		{
			int ci = vvinfo[i][j].contourID;
			int pi = vvinfo[i][j].pointID;
			points.push_back(CParticleF(contours[ci].X[pi], contours[ci].Y[pi]));
			vinfo.push_back(vvinfo[i][j]);
		}
	}

	Triangulation::Triangulator trmap(points);
	const int dims[] = {contours.size(), contours.size()};
	vector<float> W(dims[0]*dims[1], 0.0f);
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		Triangulation::_Internal::_triangle* tr = trmap.faces[i];
		for(int j=0; j<3; ++j)
		{
			Triangulation::_Internal::_edge* edge = tr->edges[j];
			int k1 = distance(points.begin(), find(points.begin(), points.end(), edge->vertices[0]->p));
			int k2 = distance(points.begin(), find(points.begin(), points.end(), edge->vertices[1]->p));
			if(vinfo[k1].contourID==vinfo[k2].contourID)
			{
				int cid = vinfo[k1].contourID;
				Triangulation::_Internal::_vertex* op = oppositeVertex(tr, edge);
				int k3 = distance(points.begin(), find(points.begin(), points.end(), op->p));
				if(vinfo[k3].contourID != cid)
				{
					float len = Distance(edge->vertices[0]->p, edge->vertices[1]->p);
					SetData2(W, vinfo[k3].contourID, cid, dims[0], dims[1], 
						GetData2(W, vinfo[k3].contourID, cid, dims[0], dims[1], 0.0f) + len);
				}
			}
		}
	}

	vector<float> F(dims[0]*dims[1], 0.0f);
	vector<int> sinks;
	for(int i=0; i<dims[1]; ++i)
	{
		bool bSink = true;
		bool bSource = true;
		for(int j=0; j<dims[0]; ++j)
		{
			float w1 = GetData2(W, j, i, dims[0], dims[1], 0.0f);
			float w2 = GetData2(W, i, j, dims[0], dims[1], 0.0f);
			SetData2(F, j, i, dims[0], dims[1], w2-w1);  //outgoing flow - incoming flow
			if(w2>w1)
			{
				bSource = false;
			}
			else if(w2<w1)
			{
				bSink = false;
			}
		}
		if(bSink)
		{
			sinks.push_back(i+1);
		}
	}

	if(nlhs >= 1)
	{
		plhs[0] = StoreData(W, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[]={sinks.size(), 1};
		plhs[2] = StoreData(sinks, mxINT32_CLASS, 2, dims);
	}

	mexUnlock();
}

