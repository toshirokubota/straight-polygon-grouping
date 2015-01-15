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
#include <Triangulation.h>
#include <FragmentInfo.h>
#include <TriangulationHelper.h>
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Kruskal.h>

GraphFactory<Triangulation::_Internal::_triangle*>* GraphFactory<Triangulation::_Internal::_triangle*>::_instance = NULL;

#include <Rank.h>
struct Treble
{
	Treble()
	{
		edge = NULL; 
		left = NULL;
		right = NULL;
		leftWeight = 0;
		rightWeight = 0;
		fitness = 0;
	}
	Triangulation::_Internal::_edge* edge;
	Triangulation::_Internal::_edge* left;
	Triangulation::_Internal::_edge* right;
	float leftWeight;
	float rightWeight;
	float fitness;
};

vector<Treble>
Saliency(Triangulation::Triangulator& trmap)
{
	vector<float> vlen(trmap.edges.size(), 0);
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		vlen[i] = trmap.edges[i]->Length();
	}
	float alpha = 10.0;
	float beta = Rank(vlen, (int)(vlen.size() * 0.05)) * 5.0;
	printf("beta = %f\n", beta);
	vector<Treble> trebles(trmap.edges.size());
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* edge = trmap.edges[i];
		float prod = 1.0-1/(1+exp(-beta*(edge->Length()-beta)/beta)); //exp(-edge->Length()*edge->Length()/(2*sigma*sigma));
		Triangulation::_Internal::_edge* chosen[] = {NULL, NULL};
		float wgts[2] = {0.0f, 0.0f};
		for(int j=0; j<2; ++j)
		{
			Triangulation::_Internal::_vertex* u = edge->vertices[j];
			Triangulation::_Internal::_vertex* v = edge->vertices[j==0 ? 1: 0]; //the other vertex
			float best = 0;
			for(int k=0; k<u->edges.size(); ++k)
			{
				Triangulation::_Internal::_edge* edge2 = u->edges[k];
				Triangulation::_Internal::_vertex* w = edge2->vertices[0]==u ? edge2->vertices[1]: edge2->vertices[0];
				if(edge == edge2) continue;
				float angle = GetVisualAngle(v->p.m_X, v->p.m_Y, w->p.m_X, w->p.m_Y, u->p.m_X, u->p.m_Y);
				float len = edge2->Length();
				float score = (1.0-cos(angle))/2.0f * (1.0-1/(1+exp(-beta*(len-beta)/beta)));
				if(score > best) 
				{
					best = score;
					chosen[j] = edge2;
					wgts[j] = (1.0-cos(angle))/2.0f;
				}
			}
			prod *= best;
		}
		trebles[i].edge = edge;
		trebles[i].left = chosen[0];
		trebles[i].right = chosen[1];
		trebles[i].fitness = prod;
		trebles[i].leftWeight = wgts[0];
		trebles[i].rightWeight = wgts[1];
	}			
	return trebles;
}

vector<CParticleF> collectCenterPoints(vector<Triangulation::_Internal::_triangle*>& faces)
{
	vector<CParticleF> center(faces.size());
	for(int i=0; i<faces.size(); ++i)
	{
		CParticleF p[3] = {faces[i]->vertices[0]->p, faces[i]->vertices[1]->p, faces[i]->vertices[2]->p};
		float x=0;
		float y=0;
		for(int j=0; j<3; ++j)
		{
			x += p[j].m_X;
			y += p[j].m_Y;
		}
		center[i] = CParticleF(x/3.0f, y/3.0f, 0.0f);
	}
	return center;
}

vector<int> 
connectedComponents(vector<Edge<Triangulation::_Internal::_triangle*>*>& forrest, 
					map<Triangulation::_Internal::_triangle*,int>& fidxmap)
{
	vector<Node<int>*> nodes;
	for(int i=0; i<fidxmap.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for(int i=0; i<forrest.size(); ++i)
	{
		int k1 = fidxmap[forrest[i]->u->key];
		int k2 = fidxmap[forrest[i]->v->key];
		merge(nodes[k1], nodes[k2]);
	}
	vector<Node<int>*> cl = clusters(nodes);
	printf("There are %d components.\n", cl.size());
	vector<int> label(nodes.size());
	for(int i=0; i<nodes.size(); ++i)
	{
		int j = distance(cl.begin(), find(cl.begin(), cl.end(), findset(nodes[i])));
		label[i] = j;
	}
	return label;
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
	float thres = 0.5;
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[1],classMode);
	} 

	Triangulation::Triangulator trmap(points);
	vector<Treble> trebles = Saliency(trmap);
	map<Triangulation::_Internal::_vertex*,int> pidxmap;
	for(int i=0; i<trmap.points.size(); ++i)
	{
		pidxmap[trmap.points[i]]=i;
	}
	map<Triangulation::_Internal::_edge*,int> eidxmap;
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		eidxmap[trmap.edges[i]]=i;
	}
	map<Triangulation::_Internal::_triangle*,int> fidxmap;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		fidxmap[trmap.faces[i]]=i;
	}
	vector<CParticleF> center = collectCenterPoints(trmap.faces);
	vector<Vertex<Triangulation::_Internal::_triangle*>*> vertices;
	vector<Edge<Triangulation::_Internal::_triangle*>*> edges;
	GraphFactory<Triangulation::_Internal::_triangle*>* factory = GraphFactory<Triangulation::_Internal::_triangle*>::GetInstance();
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(trmap.faces[i]));
	}
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		for(int j=0; j<3; ++j)
		{
			Triangulation::_Internal::_edge* edge = trmap.faces[i]->edges[j];
			Triangulation::_Internal::_triangle* tr = edge->faces[0]==trmap.faces[i] ? edge->faces[1]: edge->faces[0];
			if(tr != NULL)
			{
				int k = fidxmap[tr];
				int m = eidxmap[edge];
				Edge<Triangulation::_Internal::_triangle*>* arc = factory->makeEdge(vertices[i], vertices[k], trebles[m].fitness);
				vertices[i]->Add(arc);
				edges.push_back(arc);
			}
		}
	}
	vector<Edge<Triangulation::_Internal::_triangle*>*> mst = Kruskal(edges);
	vector<Edge<Triangulation::_Internal::_triangle*>*> forrest;
	for(int i=0; i<mst.size(); ++i)
	{
		if(mst[i]->w < thres)
		{
			forrest.push_back(mst[i]);
		}
	}
	printf("%d edges out of %d are kept with thres=%f.\n", forrest.size(), mst.size(), thres);
	vector<int> label = connectedComponents(forrest, fidxmap);
	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trmap.faces.size(), 4};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[1]->p));
			int j3 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[2]->p));
			int lb = label[fidxmap[trmap.faces[i]]];
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
			SetData2(F, i, 3, dims[0], dims[1], lb);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.edges.size(), 2};
		vector<float> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], (float)(j1+1));
			SetData2(F, i, 1, dims[0], dims[1], (float)(j2+1));
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		const int dims[] = {mst.size(), 4};
		vector<float> F(dims[0]*dims[1], 0);
		for(int i=0; i<mst.size(); ++i)
		{
			Triangulation::_Internal::_triangle* t1 = mst[i]->u->key;
			Triangulation::_Internal::_triangle* t2 = mst[i]->v->key;
			int k1 = fidxmap[t1];
			int k2 = fidxmap[t2];
			CParticleF p1 = center[k1];
			CParticleF p2 = center[k2];
			SetData2(F, i, 0, dims[0], dims[1], p1.m_X);
			SetData2(F, i, 1, dims[0], dims[1], p1.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], p2.m_X);
			SetData2(F, i, 3, dims[0], dims[1], p2.m_Y);
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 5)
	{
		const int dims[] = {forrest.size(), 4};
		vector<float> F(dims[0]*dims[1], 0);
		for(int i=0; i<forrest.size(); ++i)
		{
			Triangulation::_Internal::_triangle* t1 = forrest[i]->u->key;
			Triangulation::_Internal::_triangle* t2 = forrest[i]->v->key;
			int k1 = fidxmap[t1];
			int k2 = fidxmap[t2];
			CParticleF p1 = center[k1];
			CParticleF p2 = center[k2];
			SetData2(F, i, 0, dims[0], dims[1], p1.m_X);
			SetData2(F, i, 1, dims[0], dims[1], p1.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], p2.m_X);
			SetData2(F, i, 3, dims[0], dims[1], p2.m_Y);
		}
		plhs[4] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	factory->Clean();
	mexUnlock();
}

