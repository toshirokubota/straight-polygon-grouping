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
#include <DeformingPolygonV5.h>
using namespace StraightAxis;
int MovingParticle::_id = 0;
int DeformingPolygon::_id = 0;

struct SnapshotNode
{
	SnapshotNode(int id, int lv=0)
	{
		polygon_id = id;
		level = lv;
	}
	int polygon_id;
	int level;
	vector<CParticleF> outline;
};

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;
GraphFactory<MovingParticle*>* GraphFactory<MovingParticle*>::_instance = NULL;
GraphFactory<SnapshotNode*>* GraphFactory<SnapshotNode*>::_instance = NULL;
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
	float tthres = 100;
	if(nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(tthres,prhs[2],classMode);
	} 
	float step = tthres;
	if(nrhs >= 4)
	{
		mxClassID classMode;
		float step0;
		ReadScalar(step0,prhs[3],classMode);
		if(step0 <=1.0e-6) //too small..
		{
			printf("step = %f is too small. Reverting it to the default (%f).\n", step0, step);
		}
		else
		{
			step = step0;
		}
	} 
	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	ParticleFactory* pfactory = ParticleFactory::getInstance();

	vector<MovingParticle*> outline;
	for (int i = 0; i < points.size(); ++i)
	{
		outline.push_back(pfactory->makeParticle(points[i], Initial, 0.0f));
	}

	DeformingPolygon* dp = new DeformingPolygon(outline, 0.0f);
	vector<DeformingPolygon*> polygons = dp->deform(tthres,step);

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
		const int dims[] = {pfactory->particles.size(), 10};
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
			SetData2(F, i, 8, dims[0], dims[1], p->created);
			SetData2(F, i, 9, dims[0], dims[1], p->time);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		vector<vector<CParticleF>> contours;
		for(int i=0; i<polygons.size(); ++i)
		{
			contours.push_back(polygons[i]->snapshot);
		}
		plhs[2] = StoreContours(contours);
	}
	if(nlhs >= 4)
	{
		const int dims[] = {polygons.size(), 7};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)polygons[i]->id);
			SetData2(F, i, 1, dims[0], dims[1], (float)(polygons[i]->parent==NULL ? 0: polygons[i]->parent->id));
			SetData2(F, i, 2, dims[0], dims[1], polygons[i]->created);
			SetData2(F, i, 3, dims[0], dims[1], polygons[i]->died);
			SetData2(F, i, 4, dims[0], dims[1], (float)polygons[i]->particles.size());
			SetData2(F, i, 5, dims[0], dims[1], (float)ClockWise(polygons[i]->snapshot));
			SetData2(F, i, 6, dims[0], dims[1], (float)polygons[i]->event.type);
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 5)
	{
		vector<vector<CParticleF>> traces;
		for (int i = 0; i<polygons.size(); ++i)
		{
			vector<MovingParticle*> tr = traceBackPolygon(polygons[i]);
			vector<CParticleF> poly;
			for (int j = 0; j < tr.size(); ++j)
			{
				poly.push_back(tr[j]->p0);
			}
			traces.push_back(poly);
		}
		plhs[4] = StoreContours(traces);
	}


	//clean up
	for(int i=0; i<polygons.size(); ++i)
	{
		delete polygons[i];
	}
	GraphFactory<graphKey>::GetInstance()->Clean();
	GraphFactory<MovingParticle*>::GetInstance()->Clean();
	pfactory->clean();
	mexUnlock();
}

