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

set<pair<int,int>>
pixelize(vector<CParticleF>& polygon, int bw, int bh)
{
	set <pair<int, int>> blocks;
	if (polygon.empty()) return blocks;

	float minx, miny, maxx, maxy;
	boundingBox(polygon, minx, miny, maxx, maxy);
	int bx = (int)minx / bw;
	int by = (int)miny / bh;
	int ex = (int)maxx / bw;
	int ey = (int)maxy / bh;
	for (int i = by; i <= ey; ++i)
	{
		for (int j = bx; j <= ex; ++j)
		{
			CParticleF p((j + .5)*bw, (i + .5)*bh);
			if (inside(p, polygon))
			{
				pair<int, int> idx(j, i);
				blocks.insert(idx);
			}
		}
	}
	return blocks;
}

float 
overlapMeasure(set<pair<int, int>>& A, set<pair<int,int>>& B)
{
	set<pair<int,int>> intersect;
	set_intersection(A.begin(), A.end(), B.begin(), B.end(), std::inserter(intersect, intersect.begin()));
	set<pair<int, int>> unionset;
	set_union(A.begin(), A.end(), B.begin(), B.end(), std::inserter(unionset, unionset.begin()));
	return (float)intersect.size() / (float)unionset.size();
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
	int blockWidth = 7;
	int blockHeight = blockWidth;

	ParticleFactory* factory = ParticleFactory::getInstance();

	//convert a snapshot into a sequence of particles
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
	//compute fitness for each polygon
	vector<float> fitness(snapshots.size(), 0);
	for (int i = 0; i < polygons.size(); ++i)
	{
		fitness[i] = calculatePolygonFitness(polygons[i]);
	}
	//for each polygon, collectt a set of internal blocks
	vector<set<pair<int, int>>> pixelized(polygons.size());
	for (int i = 0; i < polygons.size(); ++i)
	{
		pixelized[i] = pixelize(polygons[i], blockWidth, blockHeight);
	}
	vector<Node<int>*> vnodes;
	for (int i = 0; i < polygons.size(); ++i)
	{
		vnodes.push_back(makeset(i));
	}
	vector<float> area(polygons.size(), 0.0f);
	for (int i = 0; i < polygons.size(); ++i)
	{
		area[i] = polygonArea(polygons[i]);
	}
	for (int i = 0; i < polygons.size(); ++i)
	{
		if (fitness[i] < thres) continue;
		for (int j = i + 1; j < snapshots.size(); ++j)
		{
			if (fitness[j] < thres) continue;
			float omeasure = overlapMeasure(pixelized[i], pixelized[j]);
			if (omeasure > overlap)
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

	if (nlhs >= 1)
	{
		const int dims[] = { reps.size(), 2 };
		vector<pair<float, int>> v;
		for (int i = 0; i < reps.size(); ++i)
		{
			int j = selected[i];
			pair<float, int> p(-fitness[j], j); //make it negative to sort in descending order
			v.push_back(p);
		}
		sort(v.begin(), v.end());
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)v[i].second+1);
			SetData2(F, i, 1, dims[0], dims[1], -v[i].first);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { fitness.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			int j = nmap[vnodes[i]->parent];
			SetData2(F, i, 0, dims[0], dims[1], fitness[i]);
			SetData2(F, i, 1, dims[0], dims[1], (float)(j + 1));
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

