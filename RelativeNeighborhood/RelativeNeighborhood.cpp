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
removeKinks(vector<Triangulation::_Internal::_vertex*>& path)
{
	vector<Triangulation::_Internal::_vertex*> res = path;
	float eps1 = 1.0e-3;
	float eps2 = 1.0e-2;
	for (int i = res.size() - 1; i >= 0 && res.size()>2; i--)
	{
		Triangulation::_Internal::_vertex* p = res[i];
		Triangulation::_Internal::_vertex* q = res[(i+1) % res.size()];
		Triangulation::_Internal::_vertex* r = res[(i-1+res.size()) % res.size()];
		float d = Distance(q->p, p->p);
		float a = GetVisualAngle(q->p.m_X, q->p.m_Y, r->p.m_X, r->p.m_Y, p->p.m_X, p->p.m_Y);
		if (d < eps1 || a < eps2)
		{
			res.erase(res.begin() + i);
			printf("removeKinks: removing a partcile.\n");
		}
	}
	return res;
}

vector<Triangulation::_Internal::_vertex*>
subsamplePath(vector<Triangulation::_Internal::_vertex*>& path, int maxGap) {
	vector<Triangulation::_Internal::_vertex*> spath;

	int nfill = (path.size() - 1) / maxGap;
	float inc = (float)(path.size()-1) / (nfill + 1);
	//printf("len=%d, gap=%d, n=%d, inc=%f\n", path.size(), maxGap, nfill, inc);
	spath.push_back(path[0]);
	float next = inc;
	for (int i = 1; i < path.size() - 1; ++i)
	{
		if (i >= next)
		{
			spath.push_back(path[i]);
			next += inc;
		}
	}
	spath.push_back(path[path.size() - 1]);
	return spath;
}

vector<Triangulation::_Internal::_vertex*>
trace2Junctions(Triangulation::_Internal::_edge* edge, Triangulation::_Internal::_vertex* source)
{
	vector<Triangulation::_Internal::_vertex*> trace;
	Triangulation::_Internal::_vertex* u = edge->vertices[0]==source ? edge->vertices[0]: edge->vertices[1];
	Triangulation::_Internal::_vertex* v = edge->vertices[0] == source ? edge->vertices[1] : edge->vertices[0];
	trace.push_back(u);
	while (v->edges.size() == 2 && v != source) 
	{
		trace.push_back(v);
		for (int i = 0; i < 2; ++i)
		{
			Triangulation::_Internal::_edge* ed = v->edges[i];
			if (ed->vertices[0] != u && ed->vertices[1] != u)
			{
				edge = ed;
				if (edge->vertices[0] == v)
				{
					v = edge->vertices[1];
					u = edge->vertices[0];
				}
				else
				{
					v = edge->vertices[0];
					u = edge->vertices[1];
				}
				break;
			}
		}
	}
	trace.push_back(v);
	return trace;
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
	float scale = 2.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(scale, prhs[1], classMode);
	}
	int maxSkip = 5; 
	if (nrhs >= 3)
	{
	mxClassID classMode;
	ReadScalar(maxSkip, prhs[2], classMode);
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
	//first, to facilitate the processing, remove non-kept edges from each vertex.
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
	for (int i = 0; i < trmap.points.size(); ++i)
	{
		Triangulation::_Internal::_vertex* p = trmap.points[i];
		if (p->edges.size()!=2)
		{
			junctions.push_back(p);
		}
	}
	//for each junction, find adjacent junctions
	set<Triangulation::_Internal::_vertex*> vkept;
	vector<pair<Triangulation::_Internal::_vertex*, Triangulation::_Internal::_vertex*>> ekept;
	for (int i = 0; i < junctions.size(); ++i)
	{
		bool bLoop = false;
		for (int j = 0; j < junctions[i]->edges.size(); ++j)
		{
			Triangulation::_Internal::_edge* ed = junctions[i]->edges[j];
			vector<Triangulation::_Internal::_vertex*> trace = trace2Junctions(ed, junctions[i]);

			//remove small traces. TK!!! - I am not sure if this is a right thing to do.
			if (trace.size() < maxSkip) continue;

			if (pmap[trace[0]] < pmap[trace[trace.size() - 1]] || (bLoop == false && pmap[trace[0]] == pmap[trace[trace.size() - 1]]))
			{
				if (pmap[trace[0]] == pmap[trace[trace.size() - 1]])
				{
					bLoop = true;
				}
				else
				{
					trace = removeKinks(subsamplePath(trace, maxSkip)); //subsample on non-loopy trace
				}

				for (int k = 0; k < trace.size(); ++k)
				{
					vkept.insert(trace[k]);
				}
				for (int k = 1; k < trace.size(); ++k)
				{
					pair<Triangulation::_Internal::_vertex*, Triangulation::_Internal::_vertex*> ed(trace[k - 1], trace[k]);
					if (ed.first != ed.second)
					{
						//because of loops, it is possible that two are the same. Do not include such degenerate cases.
						ekept.push_back(ed);
					}
				}
			}
		}
	}

	vector<Triangulation::_Internal::_vertex*> veckept;
	map<Triangulation::_Internal::_vertex*, int> mapkept;
	for (set<Triangulation::_Internal::_vertex*>::iterator it = vkept.begin(); it != vkept.end(); it++)
	{
		veckept.push_back(*it);
		mapkept[*it] = veckept.size() - 1;
	}

	if (nlhs >= 3)
	{
		const int dims[] = { veckept.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < veckept.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], veckept[i]->p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], veckept[i]->p.m_Y);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 4)
	{
		const int dims[] = { ekept.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < ekept.size(); ++i)
		{
			int k1 = mapkept[ekept[i].first];
			int k2 = mapkept[ekept[i].second];
			SetData2(F, i, 0, dims[0], dims[1], k1 + 1);
			SetData2(F, i, 1, dims[0], dims[1], k2 + 1);
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

