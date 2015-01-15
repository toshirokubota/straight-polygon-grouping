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

struct EdgePair
{
	EdgePair(Triangulation::_Internal::_edge* l=NULL, Triangulation::_Internal::_edge* r=NULL)
	{
		left = l;
		right = r;
	}
	Triangulation::_Internal::_edge* left;
	Triangulation::_Internal::_edge* right;
};

/*
find the vertex opposite to a BOUNDARY edge, e.
*/
Triangulation::_Internal::_vertex* 
oppositeSideVertex(Triangulation::_Internal::_edge* e)
{
	Triangulation::_Internal::_triangle* face = NULL;
	if(e->faces[0]!=NULL & e->faces[1]==NULL)
	{
		face = e->faces[0];
	}
	else if(e->faces[0]==NULL & e->faces[1]!=NULL)
	{
		face = e->faces[1];
	}
	else
	{
		return NULL; //this should not happen...
	}
	for(int i=0; i<3; ++i)
	{
		if(face->vertices[i] != e->vertices[0] && face->vertices[i] != e->vertices[1])
		{
			return face->vertices[i];
		}
	}
	return NULL; //this should not happen
}

vector<int> 
supportedEdges(Triangulation::Triangulator& trmap, 
			   vector<EdgePair>& edgepairs)
{
	vector<int> supported(trmap.edges.size());
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = trmap.edges[i];
		int k0 = distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), edgepairs[i].left));
		int k1 = distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), edgepairs[i].right));

		Triangulation::_Internal::_edge* e1 = edgepairs[k0].left;
		Triangulation::_Internal::_edge* e2 = edgepairs[k0].right;
		Triangulation::_Internal::_edge* e3 = edgepairs[k1].left;
		Triangulation::_Internal::_edge* e4 = edgepairs[k1].right;
		int count = 0;
		count += (e1==e ? 1: 0);
		count += (e2==e ? 1: 0);
		count += (e3==e ? 1: 0);
		count += (e4==e ? 1: 0);
		supported[i] = count;
	}
	return supported;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [P E F] = Testbed(P[, thres])");
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
	float thres = 0.8;
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[1],classMode);
	}

	Triangulation::Triangulator trmap(points);
	vector<EdgePair> edgepairs(trmap.edges.size());

	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* edge = trmap.edges[i];
		Triangulation::_Internal::_vertex* u = edge->vertices[0];
		Triangulation::_Internal::_vertex* v = edge->vertices[1];
		float mind1 = std::numeric_limits<float>::infinity();
		for(int j=0; j<u->edges.size(); ++j)
		{
			if(u->edges[j] == edge) continue;
			Triangulation::_Internal::_edge* ue = u->edges[j];
			Triangulation::_Internal::_vertex* x = ue->vertices[0]==u ? ue->vertices[1]: ue->vertices[0];
			float d = Distance(u->p, x->p);
			if(d < mind1)
			{
				mind1 = d;
				edgepairs[i].left = ue;
			}
		}
		float mind2 = std::numeric_limits<float>::infinity();
		for(int k=0; k<v->edges.size(); ++k)
		{
			if(v->edges[k] == edge) continue;
			Triangulation::_Internal::_edge* ve = v->edges[k];
			Triangulation::_Internal::_vertex* y = ve->vertices[0]==v ? ve->vertices[1]: ve->vertices[0];
			float d = Distance(v->p, y->p);
			if(d < mind2)
			{
				mind2 = d;
				edgepairs[i].right = ve;
			}
		}
	}

	vector<int> supported = supportedEdges(trmap, edgepairs);
	vector<int> saliency(trmap.edges.size(), 0);
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* edge = trmap.edges[i];
		Triangulation::_Internal::_vertex* u = edge->vertices[0];
		Triangulation::_Internal::_vertex* v = edge->vertices[1];
		if(supported[i]==1)
		{
			saliency[i] = 2; //supported from only one side
			continue;
		}
		else if(supported[i]>1)
		{
			saliency[i] = 3; //supported from both sides
			continue;
		}
		vector<EdgePair> pairs;
		for(int j=0; j<u->edges.size(); ++j)
		{
			if(u->edges[j] == edge) continue;
			Triangulation::_Internal::_edge* ue = u->edges[j];
			int jd = distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), ue));
			if(supported[jd])
			{
				for(int k=0; k<v->edges.size(); ++k)
				{
					if(v->edges[k] == edge) continue;
					Triangulation::_Internal::_edge* ve = v->edges[k];
					int kd = distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), ve));
					if(supported[kd])
					{
						pairs.push_back(EdgePair(ue, ve));
					}
				}
			}
		}
		for(int j=0; j<pairs.size(); ++j)
		{
			Triangulation::_Internal::_edge* e1 = pairs[j].left;
			Triangulation::_Internal::_edge* e2 = pairs[j].right;
			Triangulation::_Internal::_vertex* x = e1->vertices[0] == u ? e1->vertices[1]: e1->vertices[0];
			Triangulation::_Internal::_vertex* y = e2->vertices[0] == v ? e2->vertices[1]: e2->vertices[0];
			float a1 = GetVisualAngle(x->p.m_X, x->p.m_Y, v->p.m_X, v->p.m_Y, u->p.m_X, u->p.m_Y);
			float a2 = GetVisualAngle(y->p.m_X, y->p.m_Y, u->p.m_X, u->p.m_Y, v->p.m_X, v->p.m_Y);
			float sal = ((1 - cos(a1)) * (1-cos(a2)) / 4.0f);
			if(sal > thres)
			{
				saliency[i] = 1;
			}
		}
	}

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
		const int dims[] = {trmap.edges.size(), 4};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], (int)trmap.edges[i]->type+1);
			SetData2(F, i, 3, dims[0], dims[1], saliency[i]);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
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
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	mexUnlock();
}

