#include <iostream>
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
#include <szDistanceTransform.h>
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>
#include <intersectionConvexPolygons.h>
#include <szConvexHull2D.h>
#include <TriangulationSaliencyUtility.h>

float 
tensorDifference(StructuralTensor& s, StructuralTensor& t)
{
	float a1 = atan2(s.major_axis[1], s.major_axis[0]);
	float a2 = atan2(t.major_axis[1], t.major_axis[0]);
	float sn = sin(a1-a2);
	return Abs(s.major - t.major) * sn * sn;
}

vector<float>
deriveSaliency(vector<StructuralTensor>& tensors, 
				Triangulation::Triangulator& trmap)
{
	map<Triangulation::_Internal::_vertex*,int> imap;
	for(int i=0; i<trmap.points.size(); ++i)
	{
		imap[trmap.points[i]] = i;
	}

	vector<float> values(tensors.size(), 0);
	for(int i=0; i<trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* u = trmap.points[i];
		float sum = 0;
		for(int j=0; j<u->edges.size(); ++j)
		{
			Triangulation::_Internal::_vertex* v = u->edges[j]->vertices[0]==u ? u->edges[j]->vertices[1]: u->edges[j]->vertices[0];
			StructuralTensor s = tensors[imap[v]];
			sum += tensorDifference(s, tensors[i]);
		}
		sum /= u->edges.size();
		values[i] = sum;
	}
	return values;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 2 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [J] = Testbed(points, faces)");
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
	//triangle (indices to the points)
	vector<Triangulation::_Internal::_indexed_triangle> T;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
		T = indices2structs(points, T0, dimsT);
	}

	Triangulation::Triangulator trmap(points, T);
	vector<StructuralTensor> tensors(trmap.points.size());
	for(int i=0; i<trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* u = trmap.points[i];
		vector<CParticleF> pnts(u->edges.size());
		for(int j=0; j<u->edges.size(); ++j)
		{
			Triangulation::_Internal::_vertex* v = u->edges[j]->vertices[0]==u ? u->edges[j]->vertices[1]: u->edges[j]->vertices[0];
			float len = u->edges[j]->Length();
			pnts[j].m_X = (v->p.m_X - u->p.m_X); // / (len*len);
			pnts[j].m_Y = (v->p.m_Y - u->p.m_Y); // / (len*len);
		}
		tensors[i] = StructuralTensor(pnts);
		if(i+1==198 || i+1==228 || i+1==28)
		{
			printf("%d::\n", i+1);
			for(int j=0; j<pnts.size(); ++j)
			{
				printf("%f %f\n", pnts[j].m_X, pnts[j].m_Y);
			}
			printf("%f %f %f %f %f %f %f %f\n", 
				tensors[i].cx, tensors[i].cy, tensors[i].major, tensors[i].minor, 
				tensors[i].major_axis[0], tensors[i].major_axis[1], tensors[i].minor_axis[0], tensors[i].minor_axis[1]);
			//break;
		}
		tensors[i].cx += u->p.m_X;
		tensors[i].cy += u->p.m_Y;
	}
	vector<float> saliency = deriveSaliency(tensors, trmap);
	
	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 8};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], tensors[i].cx);
			SetData2(F, i, 1, dims[0], dims[1], tensors[i].cy);
			SetData2(F, i, 2, dims[0], dims[1], tensors[i].major);
			SetData2(F, i, 3, dims[0], dims[1], tensors[i].minor);
			SetData2(F, i, 4, dims[0], dims[1], tensors[i].major_axis[0]);
			SetData2(F, i, 5, dims[0], dims[1], tensors[i].major_axis[1]);
			SetData2(F, i, 6, dims[0], dims[1], tensors[i].minor_axis[0]);
			SetData2(F, i, 7, dims[0], dims[1], tensors[i].minor_axis[1]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {saliency.size(), 1};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], saliency[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

