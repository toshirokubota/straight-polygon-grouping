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
#include <MovingParticle.h>
#include <ParticleSimulator.h>

GraphFactory<CParticleF>* GraphFactory<CParticleF>::_instance = NULL;
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
	ParticleSimulator simulator;
	//Points
	vector<CParticleF> points; 
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		if (dimsP[1] == ParticleDumpSize)
		{
			float time0 = std::numeric_limits<float>::infinity();
			simulator.LoadParticles(P0, dimsP);
			if (nrhs >= 2)
			{
				mxClassID classMode;
				vector<Snapshot> snapshots = Snapshot::LoadSnapshots(prhs[1]);
				simulator.Restore(snapshots);
			}
		}
		else
		{
			for (int i = 0; i < dimsP[0]; ++i)
			{
				float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
				float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
				points.push_back(CParticleF(x, y));
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
			simulator.Prepare(points, E); 
		}
	}
	float endtime = 2.0f;
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
	bool bdebug = false;
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(bdebug, prhs[4], classMode);
	}

	simulator.Simulate(endtime, delta, bdebug); 

	if (nlhs >= 1) 
	{
		plhs[0] = Snapshot::StoreSnapshots(simulator.snapshots);
	}
	if (nlhs >= 2)
	{
		plhs[1] = simulator.SaveParticles();
	}
	if (nlhs >= 3)
	{
		plhs[2] = Snapshot::StoreSnapshots(simulator.closedRegions);
	}
	if (nlhs >= 4)
	{
		plhs[3] = simulator.SaveDoneEvents();
	}

	ParticleFactory::getInstance()->clean();
	mexUnlock();
}

