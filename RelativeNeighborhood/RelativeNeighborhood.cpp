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
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>

vector<Triangulation::_Internal::_vertex*> 
trace2Junctions(Triangulation::_Internal::_vertex* p,
set<Triangulation::_Internal::_vertex*>& junctions)
{
	vector<Triangulation::_Internal::_vertex*> neighbors;
	//do BFS until it hits a junction
	set<Triangulation::_Internal::_vertex*> pset;
	pset.insert(p);
	vector<Triangulation::_Internal::_vertex*> Q(1, p);
	while (Q.empty() == false)
	{
		vector<Triangulation::_Internal::_vertex*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			Triangulation::_Internal::_vertex* q = Q[i];
			for (int j = 0; j < q->edges.size(); ++j)
			{
				Triangulation::_Internal::_edge* e = q->edges[j];
				Triangulation::_Internal::_vertex* r = e->vertices[0] == q ? e->vertices[1] : e->vertices[0];
				if (pset.find(r) == pset.end())
				{
					if (junctions.find(r) != junctions.end())
					{
						neighbors.push_back(r);
					}
					else
					{
						Q2.push_back(r);
					}
					pset.insert(r);
				}
			}
		}
		Q = Q2;
	}
	return neighbors;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "TriangleTree", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [P F E T] = TriangleTree(P)");
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
	float thres = 10.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
	}
	float scale = 2.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(scale, prhs[2], classMode);
	}

	Triangulation::Triangulator trmap(points);

	map<Triangulation::_Internal::_vertex*,int> pmap;
	for(int i=0; i<trmap.points.size(); ++i)
	{
		pmap[trmap.points[i]] = i;
	}
	map<Triangulation::_Internal::_edge*,int> emap;
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		emap[trmap.edges[i]] = i;
	}
	vector<float> W(trmap.edges.size());
	for (int i = 0; i<trmap.edges.size(); ++i)
	{
		W[i] = trmap.edges[i]->Length();
	}

	//compute the relative neighborhood graph
	vector<bool> keep(trmap.edges.size(), true);
	for (int i = 0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = trmap.edges[i];
		Triangulation::_Internal::_vertex* u = e->vertices[0];
		Triangulation::_Internal::_vertex* v = e->vertices[1];
		if (W[i]>thres)
		{
			keep[i] = false;
			continue;;
		}
		for (int k1 = 0; k1 < u->edges.size() && keep[i]; ++k1)
		{
			Triangulation::_Internal::_edge* e2 = u->edges[k1];
			Triangulation::_Internal::_vertex* w = e2->vertices[0] == u ? e2->vertices[1] : e2->vertices[0];
			for (int k2 = 0; k2 < v->edges.size() && keep[i]; ++k2)
			{
				Triangulation::_Internal::_edge* e3 = v->edges[k2];
				Triangulation::_Internal::_vertex* x = e3->vertices[0] == v ? e3->vertices[1] : e3->vertices[0];
				if (w == x)
				{
					if (W[i] > max(W[emap[e2]], W[emap[e3]]))
					{
						keep[i] = false;
					}
				}
			}
		}
	}
	//further remove those that are too long compared to two least ones
	//first get non-essential edges for each vertex, which are candidate for removal
	vector<vector<Triangulation::_Internal::_edge*>> ness(trmap.points.size());
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = trmap.points[i];
		vector<float> vlen;
		for (int j = 0; j < p->edges.size(); ++j) {
			Triangulation::_Internal::_edge* e = p->edges[j];
			int k = emap[e];
			if (keep[k] == true)
			{
				vlen.push_back(W[k]);
			}
		}
		if (vlen.size() < 2) //continue; //this should not happen, I think.
		{
			mexWarnMsgTxt("RelativeNeighborhood: dangling edge found.");
			continue;
		}

		sort(vlen.begin(), vlen.end());
		float thres = vlen[1] * scale;
		vector<Triangulation::_Internal::_edge*> ne;
		int count=0;
		for (int j = 0; j < p->edges.size(); ++j) {
			Triangulation::_Internal::_edge* e = p->edges[j];
			int k = emap[e];
			if (keep[k] == true)
			{
				count++;
				if (W[k] > thres)
				{
					ne.push_back(e);
				}
			}
		}
		ness[i] = ne;
	}
	//now remove any kept edges that are deemed non-essential by both end vertices.
	for (int i = 0; i<trmap.edges.size(); ++i)
	{
		if (keep[i] == false) continue;
		Triangulation::_Internal::_edge* e = trmap.edges[i];
		Triangulation::_Internal::_vertex* u = e->vertices[0];
		Triangulation::_Internal::_vertex* v = e->vertices[1];
		int ui = pmap[u];
		int vi = pmap[v];
		if (find(ness[ui].begin(), ness[ui].end(), e) != ness[ui].end() &&
			find(ness[vi].begin(), ness[vi].end(), e) != ness[vi].end())
		{
			keep[i] = false;
		}
	}

	if (nlhs >= 1)
	{
		const int dims[] = { trmap.points.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { trmap.edges.size(), 3 };
		vector<int> F(dims[0] * dims[1], 0);
		for (int i = 0; i<trmap.edges.size(); ++i)
		{
			int j1 = pmap[trmap.edges[i]->vertices[0]];
			int j2 = pmap[trmap.edges[i]->vertices[1]];
			SetData2(F, i, 0, dims[0], dims[1], (j1 + 1));
			SetData2(F, i, 1, dims[0], dims[1], (j2 + 1));
			SetData2(F, i, 2, dims[0], dims[1], keep[i] ? 1 : 0);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	//finally reduce the graph by eliminating vertices that are not junctions (ones with more than 2 edges)
	//first facilitate the processing, remove non-kept edges from each vertex.
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = trmap.points[i];
		for (int j = p->edges.size() - 1; j >= 0; j--)
		{
			int k = emap[p->edges[j]];
			if (keep[k] == false)
			{
				p->edges.erase(p->edges.begin() + j);
			}
		}
	}
	vector<Triangulation::_Internal::_vertex*> junctions;
	set<Triangulation::_Internal::_vertex*> jset;
	map<Triangulation::_Internal::_vertex*,int> jmap;
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = trmap.points[i];
		if (p->edges.size()>2)
		{
			junctions.push_back(p);
			jset.insert(p);
			jmap[p] = junctions.size() - 1;
		}
	}
	//for each junction, find adjacent junctions
	vector<pair<int, int>> newlinks; //collect all edges without duplicates
	for (int i = 0; i < junctions.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = junctions[i];
		vector<Triangulation::_Internal::_vertex*> neighbors = trace2Junctions(p, jset);
		int k1 = jmap[p];
		for (int j = 0; j < neighbors.size(); ++j)
		{
			int k2 = jmap[neighbors[j]];
			if (k1 < k2)
			{
				newlinks.push_back(pair<int, int>(k1, k2));
			}
		}
	}

	if (nlhs >= 3)
	{
		const int dims[] = { junctions.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], junctions[i]->p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], junctions[i]->p.m_Y);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { newlinks.size(), 2 };
		vector<int> F(dims[0] * dims[1], 0);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], newlinks[i].first+1);
			SetData2(F, i, 1, dims[0], dims[1], newlinks[i].second + 1);
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

