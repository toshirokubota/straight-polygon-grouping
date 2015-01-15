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
#include <hash_map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <DisjointSet.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
#include <StraightOffsetPolygon.h>

typedef CParticleF graphKey;

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;
GraphFactory<MovingParticle*>* GraphFactory<MovingParticle*>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;

int MovingParticle::_id = 0;

mxArray*
StoreContours(const vector<vector<CParticleF>>& polygons)
{
	const int dims[] = {polygons.size()};
	mxArray* cell = mxCreateCellArray(1, (mwSize*) dims);
	for(int i=0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size() ;
		const int dimsC[] = {n, 3};
		mxArray* ar = mxCreateNumericArray(2, (mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = polygons[i][j].m_X;
			p[n+j] = polygons[i][j].m_Y;
			p[2*n+j] = polygons[i][j].m_Z;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

vector<pair<int,int>>
indices2pairs(vector<int> T, const int* dims)
{
	int n = dims[0];
	vector<pair<int,int>> pairs(n);
	for(int i=0; i<n; ++i)
	{
		int k1 = GetData2(T, i, 0, dims[0], dims[1], 0);
		int k2 = GetData2(T, i, 1, dims[0], dims[1], 0);
		pairs[i] = pair<int,int>(k1-1, k2-1);
	}
	return pairs;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	//Points
	vector<CParticleF> points;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(CParticleF(x, y));
			//points.insert(points.begin(), CParticleF(x,y));
		}
	}
	//edges (indices to the points)
	vector<pair<int, int>> E;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
		E = indices2pairs(T0, dimsT);
	}
	float step = 1.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(step, prhs[2], classMode);
	}
	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	ParticleFactory* pfactory = ParticleFactory::getInstance();

	vector<Vertex<graphKey>*> vertices;
	for (int i = 0; i < points.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(points[i]));
	}
	vector<Edge<graphKey>*> edges;
	for (int i = 0; i < E.size(); ++i)
	{
		pair<int, int> idx = E[i];
		Edge<graphKey>* edge = factory->makeEdge(vertices[idx.first], vertices[idx.second], 1.0f);
		Edge<graphKey>* edge2 = factory->makeEdge(vertices[idx.second], vertices[idx.first], 1.0f);
		vertices[idx.first]->Add(edge);
		vertices[idx.second]->Add(edge2);
		edges.push_back(edge);
		edges.push_back(edge2);
	}
	vector<StraightOffsetPolygon*> polygons = traceForrest(vertices);
	float time = 0.05f;
	for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
	{
		(*it)->update(time);
	}
	for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
	{
		(*it)->updateEvent();
	}
	for (int iter = 0; iter < 10; ++iter)
	{
		MovingParticle* particle = pfactory->getNext();
		if (particle == NULL) break;
		particle->applyEvent();
		float delta = particle->getEvent().t - time;
		for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
		{
			(*it)->update(delta);
		}
		for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
		{
			(*it)->updateEvent();
		}
		time = particle->getEvent().t;
	}
	if(nlhs >= 1)
	{
		vector<vector<float>> F;
		const int dims[] = {polygons.size(), 1};
		for (int i = 0; i < polygons.size(); ++i)
		{
			const int dims2[] = { polygons[i]->size(), 2 };
			vector<float> P(dims2[0]*dims2[1]);
			for (int j = 0; j < polygons[i]->size(); ++j)
			{
				MovingParticle* p = polygons[i]->get(j);
				SetData2(P, j, 0, dims2[0], dims2[1], p->getP().m_X);
				SetData2(P, j, 1, dims2[0], dims2[1], p->getP().m_Y);
			}
			F.push_back(P);
		}
		plhs[0] = StoreDataCell(F, mxSINGLE_CLASS, 2, dims, 2);
	}

	for (int i = 0; i < polygons.size(); ++i)
	{
		delete polygons[i];
	}

	GraphFactory<graphKey>::GetInstance()->Clean();
	GraphFactory<MovingParticle*>::GetInstance()->Clean();
	pfactory->clean();
	mexUnlock();
}

