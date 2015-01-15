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
#include <FragmentInfo.h>
#include <DotGroupUtility.h>
#include <TriangulationHelper.h>
#include <MovingParticle.h>
#include <ParticleSimulator.h>
#include <Snapshot.h>
#include <PolygonAreaUtility.h>
#include <Graph.h>
#include <GraphFactory.h>

GraphFactory<CParticleF>* GraphFactory<CParticleF>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;
int MovingParticle::_id = 0;

float
calculatePolygonFitness(vector<CParticleF>& polygon)
{
	float area = polygonArea(polygon);
	float len2 = 0;
	for (int i = 0; i < polygon.size(); ++i)
	{
		int j = i == polygon.size()-1 ? 0 : i + 1;
		float d = Distance(polygon[i], polygon[j]);
		len2 += d * d;
	}
	return area / len2;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: areas = OverlapArea(P, I, Q, J)");
		return;
	}
	vector<float> P0;
	const int* dimsP;
	mxClassID classIdP;
	int ndimP;
	LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
	ParticleSimulator simulator;
	simulator.LoadParticles(P0, dimsP);

	vector<Snapshot> snapshots = Snapshot::LoadSnapshots(prhs[1]);
	float thres = 0.5f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[2], classMode);
	}
	float overlap = 0.5f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(overlap, prhs[3], classMode);
	}

	ParticleFactory* factory = ParticleFactory::getInstance();

	vector<vector<CParticleF>> polygons(snapshots.size());
	for (int i = 0; i < snapshots.size(); ++i)
	{
		vector<CParticleF> polygon;
		for (int j = 0; j < snapshots[i].size(); ++j)
		{
			int id = snapshots[i].get(j);
			MovingParticle* p = factory->get(id);
			polygon.push_back(p->project(snapshots[i].getTime()));
		}
		polygons[i] = polygon;
	}

	vector<float> fitness(snapshots.size(), 0);
	for (int i = 0; i < polygons.size(); ++i)
	{
		fitness[i] = calculatePolygonFitness(polygons[i]);
	}

	vector<Node<int>*> vnodes;
	for (int i = 0; i < polygons.size(); ++i)
	{
		vnodes.push_back(makeset(i));
	}
	vector<vector<vector<CParticleF>>> triangulated(polygons.size()); //oh boy. for the first time vector of vector of vector...
	vector<float> area(polygons.size(), 0.0f);
	for (int i = 0; i < polygons.size(); ++i)
	{
		bool b;
		triangulated[i] = triangulateInsidePolygon(polygons[i], b);
		if (b)
		{
			area[i] = polygonArea(triangulated[i]);
		}
	}
	for (int i = 0; i < polygons.size(); ++i)
	{
		if (fitness[i] < thres) continue;
		for (int j = i + 1; j < snapshots.size(); ++j)
		{
			if (fitness[j] < thres) continue;
			vector<vector<CParticleF>> overlap = polygonOverlap(triangulated[i], triangulated[j]);
			float area_o = polygonArea(overlap);
			if (area_o / (area[i] + area[j]) > thres)
			{
				merge(vnodes[i], vnodes[j]);
			}
		}
	}

	vector<Node<int>*> reps = clusters(vnodes);

	vector<int> selected(reps.size());
	map<Node<int>*, int> nmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		selected[i] = reps[i]->key;
		nmap[reps[i]] = i;
	}
	for (int i = 0; i < vnodes.size(); ++i)
	{
		int j = nmap[vnodes[i]->parent];
		if (fitness[i] > fitness[selected[j]])
		{
			selected[j] = i;
		}
	}

	if(nlhs >= 1)
	{
		const int dims[] = { fitness.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			int j = nmap[vnodes[i]->parent];
			SetData2(F, i, 0, dims[0], dims[1], fitness[i]);
			SetData2(F, i, 1, dims[0], dims[1], (float)(j+1));
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

