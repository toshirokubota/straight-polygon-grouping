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
//#include <StraightOffsetPolygon.h>
#include <MovingParticle.h>


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

bool
checkForCollision(MovingParticle* p)
{
	CParticleF p0 = p->getP0();
	CParticleF p1 = p->getP();
	CParticleF q0 = p->getNext()->getP0();
	CParticleF q1 = p->getNext()->getP();
	pair<float, float> param = _IntersectConvexPolygon::intersect(p0, p1, q0, q1);
	if (param.first > 0 && param.first <= 1.0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

MovingParticle* 
checkForSplit(MovingParticle* p)
{
	CParticleF p0 = p->getP0();
	CParticleF p1 = p->getP();
	ParticleFactory* factory = ParticleFactory::getInstance();
	set<MovingParticle*>::iterator it = factory->activeSet.begin();
	float t = numeric_limits<float>::infinity();
	MovingParticle* z = NULL;
	for (; it != factory->activeSet.end(); ++it)
	{
		MovingParticle* r = *it;
		CParticleF q0 = r->getP();
		CParticleF q1 = r->getNext()->getP();
		pair<float, float> param = _IntersectConvexPolygon::intersect(p0, p1, q0, q1);
		if (param.first > 0 && param.first <= 1.0)
		{
			if (param.first < t)
			{
				t = param.first;
				z = r;
			}
		}
	}
	return z;
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
	float endtime = 2.0f;
	//float step = 1.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(endtime, prhs[2], classMode);
	}
	float delta = 0.1f;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(delta, prhs[3], classMode);
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
	
	vector<vector<CParticleF>> snapshots;
	MovingParticle::traceForrest(vertices);
	float time = 0.0f;
	float delta0 = 0.01f;
	for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
	{
		(*it)->update(delta0);
	}
	time += delta0;
	float snapTime = delta;
	vector<EventStruct> doneEvents;
	while (time < endtime)
	{
		vector<EventStruct> events;
		for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
		{
			(*it)->updateEvent();
			events.push_back((*it)->getEvent());
		}
		sort(events.begin(), events.end());

		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;
		if (p->getId() == 16)
		{
			time += 0;
		}
		if (p->getEvent().t > endtime) break;

		for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
		{
			(*it)->update(p->getEvent().t - time);
		}
		time = p->getEvent().t;

		p->getEvent().print();
		if (p->applyEvent() == false) break;

		MovingParticle::removeUnstable();
		if (time >= snapTime)
		{
			vector<vector<CParticleF>> shots = MovingParticle::takeSnapshots();
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		MovingParticle::quickFinish();
		/*if (MovingParticle::sanityCheck() == false)
		{
			printf("Violation of sanity check found at %f.\n", time);
			vector<vector<CParticleF>> shots = MovingParticle::takeSnapshots();
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			break;
		}*/
		doneEvents.push_back(p->getEvent());
	}
	if(nlhs >= 1)
	{
		vector<vector<float>> F;
		const int dims[] = { snapshots.size(), 1 };
		for (int i = 0; i < snapshots.size(); ++i)
		{
			const int dims2[] = { snapshots[i].size(), 2 };
			vector<float> P(dims2[0]*dims2[1]);
			for (int j = 0; j < snapshots[i].size(); ++j)
			{
				CParticleF p = snapshots[i][j];
				SetData2(P, j, 0, dims2[0], dims2[1], p.m_X);
				SetData2(P, j, 1, dims2[0], dims2[1], p.m_Y);
			}
			F.push_back(P);
		}
		plhs[0] = StoreDataCell(F, mxSINGLE_CLASS, 2, dims, 2);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { pfactory->particles.size(), 10 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			MovingParticle* p = pfactory->particles[i];
			SetData2(F, i, 0, dims[0], dims[1], p->getP0().m_X);
			SetData2(F, i, 1, dims[0], dims[1], p->getP0().m_Y);
			SetData2(F, i, 2, dims[0], dims[1], p->getP().m_X);
			SetData2(F, i, 3, dims[0], dims[1], p->getP().m_Y);
			float vx, vy;
			p->getVelocity(vx, vy);
			SetData2(F, i, 4, dims[0], dims[1], vx);
			SetData2(F, i, 5, dims[0], dims[1], vy);
			SetData2(F, i, 6, dims[0], dims[1], (float)p->getId());
			SetData2(F, i, 7, dims[0], dims[1], (float)p->getType());
			SetData2(F, i, 8, dims[0], dims[1], p->getCreatedTime());
			SetData2(F, i, 9, dims[0], dims[1], p->getTime());
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	GraphFactory<graphKey>::GetInstance()->Clean();
	GraphFactory<MovingParticle*>::GetInstance()->Clean();
	pfactory->clean();
	mexUnlock();
}

