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
			//points.push_back(CParticleF(x, y));
			points.insert(points.begin(), CParticleF(x,y));
		}
	}
	int K = 0;
	if(nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(K,prhs[1],classMode);
	} 

	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	ParticleFactory* pfactory = ParticleFactory::getInstance();

	vector<DeformingPolygon> polygons;
	polygons.push_back(DeformingPolygon(points));
	int count = 0;
	bool bDone = false;
	vector<vector<CParticleF>> contours;
	while(!bDone && polygons.empty() == false)
	{
		DeformingPolygon poly = polygons[0];
		polygons.erase(polygons.begin());
		while(true)
		{
			int ec, el;
			count++;
			if(K == count)
			{
				bDone = true;
				break;
			}
			vector<DeformingPolygon> newpolys = fixPolygon(poly);
			if(poly.sanityCheck(ec, el) == false)
			{
				ec += 0;
			}
			polygons.insert(polygons.end(), newpolys.begin(), newpolys.end());
			if(poly.particles.size()>2)
			{
				float t1 = poly.nextEdgeEvent();
				float t2 = poly.nextSplitEvent();
				float t = Min(t1, t2);
				if(t >= std::numeric_limits<float>::infinity()) break;
				if(t1 > t2) printf("Split event at %d.\n", count);
				else printf("Edge event at %d.\n", count);

				poly.deform(t);
				vector<CParticleF> pnts;
				for(int k=0; k<poly.particles.size(); ++k)
				{
					pnts.push_back(poly.particles[k]->p);
				}
				contours.push_back(pnts);
			}
			else break;
			if(poly.sanityCheck(ec, el) == false)
			{
				ec += 0;
			}
		}
	}
	//vector<DeformingPolygon> newpolys = fixPolygon(polygons[0]);
	//polygons[0].deform(5.0);

	if(nlhs >= 1)
	{
		const int dims[] = {pfactory->particles.size(), 7};
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
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StorePolygons(contours);
	}
	factory->Clean();
	pfactory->clean();
	mexUnlock();
}

