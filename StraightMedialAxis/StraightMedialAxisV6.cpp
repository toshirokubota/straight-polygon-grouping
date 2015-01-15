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
typedef CParticleF graphKey;
#include <DeformingPolygonV3.h>
using namespace StraightAxis;
int MovingParticle::_id = 0;
int DeformingPolygon::_id = 0;

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;
GraphFactory<MovingParticle*>* GraphFactory<MovingParticle*>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;

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

void
printTree(MovingParticle* particle, 
          set<MovingParticle*>& traced,
		  string tab)
{
	if(traced.find(particle) == traced.end())
	{
		traced.insert(particle);
		printf("%s%d(%2.2f, %2.2f)-(%2.2f, %2.2f) %d\n", 
			tab.c_str(), particle->id, particle->p0.m_X, particle->p0.m_Y, particle->p.m_X, particle->p.m_Y, particle->type);
		for(int i=0; i<particle->descendents.size(); ++i)
		{
			printTree(particle->descendents[i], traced, tab + "\t");
		}
	}
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
		for(int i=0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(CParticleF(x, y));
			//points.insert(points.begin(), CParticleF(x,y));
		}
	}
	//triangle (indices to the points)
	vector<pair<int,int>> edges;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
		edges = indices2pairs(T0, dimsT);
	}
	float tthres = 10;
	if(nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(tthres,prhs[2],classMode);
	} 
	float delta = 0.1; 
	if(nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(delta,prhs[3],classMode);
	} 
	float eps = 0.01;
	if(nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(eps,prhs[4],classMode);
		if(eps <=1.0e-6) //too small..
		{
			printf("eps = %f is too small. Reverting it to the default (0.01).\n", eps);
			eps = 0.01;
		}
	} 
	float offset = 0.4f;
	if(nrhs >= 6)
	{
		mxClassID classMode;
		ReadScalar(offset,prhs[5],classMode);
	}
	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	ParticleFactory* pfactory = ParticleFactory::getInstance();

	vector<Vertex<graphKey>*> vertices;
	for(int i=0; i<points.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(points[i]));
	}
	for(int i=0; i<edges.size(); ++i)
	{
		pair<int,int> idx = edges[i];
		//printf("Edge %d - %d\n", idx.first+1, idx.second+1);
		Edge<graphKey>* edge = factory->makeEdge(vertices[idx.first], vertices[idx.second], 1.0f);
		Edge<graphKey>* edge2 = factory->makeEdge(vertices[idx.second], vertices[idx.first], 1.0f);
		vertices[idx.first]->Add(edge);
		vertices[idx.second]->Add(edge2);
	}
	vector<MovingParticle*> trace = traceTree(vertices);
	vector<MovingParticle*> outline = offsetPolygon(trace, offset, 1.0e-5);

	DeformingPolygon dp(0.0f);
	dp.particles = outline;
	vector<DeformingPolygon> polygons;
	polygons.push_back(dp);
	
	vector<vector<CParticleF>> contours;
	float time = 0;
	vector<float> vtimes(1, 0.0f);
	int count = 0;
	while(polygons.empty() == false && time < tthres)
	{
		count++;
		float mint = delta; //the largest amount that can be moved.
		string eventType = "Unknown";
		for(int i=0; i<polygons.size(); ++i)
		{
			pair<float,MovingParticle*> e1 = polygons[i].nextEdgeEvent();
			pair<float,MovingParticle*> e2 = polygons[i].nextSplitEvent();
			if(mint > e1.first)
			{
				mint = e1.first;
				eventType = "Edge";
			}
			if(mint > e2.first)
			{
				mint = e2.first;
				eventType = "Split";
			}
		}
		time += mint;
		printf("%d: %s %f %f %d\n", count, eventType.c_str(), mint, time, polygons.size());
		vtimes.push_back(time);
		vector<DeformingPolygon> newpolygons;
		for(int i=0; i<polygons.size(); ++i)
		{
			polygons[i].deform(mint);
			vector<DeformingPolygon> polygons2 = fixPolygon(polygons[i], time, eps, false);
			newpolygons.insert(newpolygons.end(), polygons2.begin(), polygons2.end());
		}

		vector<CParticleF> pnts;
		for(int i=0; i<newpolygons.size(); ++i)
		{
			for(int k=0; k<newpolygons[i].particles.size(); ++k)
			{
				CParticleF p = newpolygons[i].particles[k]->p;
				p.m_Z = (float)i; // use Z to separate polygons
				pnts.push_back(p);
			}
		}
		contours.push_back(pnts);

		for(int i=newpolygons.size()-1; i>=0; i--)
		{
			if(newpolygons[i].particles.size() <= 3)
			{
				quickFinish(newpolygons[i].particles, time);
				newpolygons.erase(newpolygons.begin()+i);
			}
		}
		polygons = newpolygons;
	}

	if(nlhs >= 1)
	{
		const int dims[] = {outline.size(), 4};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], outline[i]->p0.m_X);
			SetData2(F, i, 1, dims[0], dims[1], outline[i]->p0.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], outline[i]->v[0]);
			SetData2(F, i, 3, dims[0], dims[1], outline[i]->v[1]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trace.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], trace[i]->p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], trace[i]->p.m_Y);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {pfactory->particles.size(), 8};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			MovingParticle* p = pfactory->particles[i];
			SetData2(F, i, 0, dims[0], dims[1], p->p0.m_X);
			SetData2(F, i, 1, dims[0], dims[1], p->p0.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], p->p.m_X);
			SetData2(F, i, 3, dims[0], dims[1], p->p.m_Y);
			SetData2(F, i, 4, dims[0], dims[1], p->v[0]);
			SetData2(F, i, 5, dims[0], dims[1], p->v[1]);
			SetData2(F, i, 6, dims[0], dims[1], (float)p->id);
			SetData2(F, i, 7, dims[0], dims[1], (float)p->type);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreContours(contours);
	}
	if(nlhs >= 5)
	{
		plhs[4] = StoreStraightAxes(polygons);
	}

	//clean up
	GraphFactory<graphKey>::GetInstance()->Clean();
	GraphFactory<MovingParticle*>::GetInstance()->Clean();
	pfactory->clean();
	mexUnlock();
}

