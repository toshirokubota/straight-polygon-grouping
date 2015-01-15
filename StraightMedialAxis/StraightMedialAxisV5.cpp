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
#include <TriangulationWithMemento.h>
#include <FragmentInfo.h>
//#include <TriangulationHelper.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
typedef CParticleF graphKey;
#include <DeformingPolygon.h>
using namespace StraightAxis;

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;
int MovingParticle::_id = 0;
ParticleFactory* ParticleFactory::_instance = NULL;

mxArray*
StorePolygons(const vector<vector<CParticleF>>& polygons)
{
	const int dims[] = {polygons.size()};
	mxArray* cell = mxCreateCellArray(1, (mwSize*) dims);
	for(int i=0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size() ;
		const int dimsC[] = {n, 2};
		mxArray* ar = mxCreateNumericArray(2, (mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = polygons[i][j].m_X;
			p[n+j] = polygons[i][j].m_Y;
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
	float tthres = 100;
	if(nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(tthres,prhs[2],classMode);
	} 
	int K=1000; //large enough?
	if(nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(K,prhs[3],classMode);
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
	vector<CParticleF> trace = traceTree(vertices);
	vector<CParticleF> outline = offsetPolygon(trace, 0.4f, 1.0e-5);

	DeformingPolygon dp(outline);
	vector<DeformingPolygon> polygons;
	polygons.push_back(dp);
	
	int count = 0;
	int poly_count = 1;
	vector<vector<CParticleF>> contours;
	contours.push_back(outline);
	int poly_id = 1;
	while(polygons.empty() == false)
	{
		DeformingPolygon poly = polygons[0];
		polygons.erase(polygons.begin());
		float tsum = 0;
		printf("Polygon %d.\n", poly_id++);
		while(tsum < tthres && count < K)
		{
			count++;

			vector<DeformingPolygon> newpolys = fixPolygon(poly, eps, false);
			int ec = poly.sanityCheck();
			if(ec != 0)
			{
				printf("Kept:\n");
				for(int i=0; i<poly.particles.size(); ++i)
				{
					MovingParticle* p = poly.particles[i];
					printf("%d %f %f %d (%x %x %x) (%d, %d, %d)\n", 
						i+1, p->p.m_X, p->p.m_Y, p->id, p->prev, p, p->next, p->prev->id, p->id, p->next->id); 
				}
				for(int j=0; j<newpolys.size(); ++j)
				{
					printf("new poly %d\n", j+1);
					for(int i=0; i<newpolys[j].particles.size(); ++i)
					{
						MovingParticle* p = newpolys[j].particles[i];
						printf("%d %f %f %d (%x %x %x)\n", i+1, p->p.m_X, p->p.m_Y, p->id, p->prev, p, p->next); 
					}
				}
				mexErrMsgTxt("sanity checked failed.\n");
			}
			poly_count += newpolys.size();

			polygons.insert(polygons.end(), newpolys.begin(), newpolys.end());
			if(poly.particles.size() <= 3)
			{
				quickFinish(poly.particles);
				break;
			}
			else
			{
				pair<float,MovingParticle*> e1 = poly.nextEdgeEvent();
				pair<float,MovingParticle*> e2 = poly.nextSplitEvent();
				float t = Min(e1.first, e2.first);
				if(t >= std::numeric_limits<float>::infinity()) break;
				poly.deform(t);

				if(e1.first > e2.first) 
					printf("Split event at %d, %3.3f by %d (%3.3f,%3.3f)->(%3.3f,%3.3f) %f.\n", 
						count, tsum, e2.second->id+1, e2.second->p0.m_X, e2.second->p0.m_Y, e2.second->p.m_X, e2.second->p.m_Y, e2.first);
				else 
					printf("Edge event at %d, %3.3f by %d (%3.3f,%3.3f)->(%3.3f,%3.3f) %3f.\n", 
						count, tsum, e1.second->id+1, e1.second->p0.m_X, e1.second->p0.m_Y, e1.second->p.m_X, e1.second->p.m_Y, e1.first);

				vector<CParticleF> pnts;
				for(int k=0; k<poly.particles.size(); ++k)
				{
					pnts.push_back(poly.particles[k]->p);
				}
				contours.push_back(pnts);
				tsum += t;
			}
		}
	}

	if(nlhs >= 1)
	{
		const int dims[] = {outline.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], outline[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], outline[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trace.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], trace[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], trace[i].m_Y);
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
		plhs[3] = StorePolygons(contours);
	}
	factory->Clean();
	pfactory->clean();
	mexUnlock();
}

