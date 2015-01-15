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

GraphFactory<Triangulation::_Internal::_edge*>* GraphFactory<Triangulation::_Internal::_edge*>::_instance = NULL;

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

vector<int> 
connectedComponents(vector<Vertex<Triangulation::_Internal::_edge*>*>& vertices, 
					vector<Edge<Triangulation::_Internal::_edge*>*>& forrest) 
{
	map<Triangulation::_Internal::_edge*,int> idxmap;
	vector<Node<Triangulation::_Internal::_edge*>*> nodes;
	for(int i=0; i<vertices.size(); ++i)
	{
		nodes.push_back(makeset(vertices[i]->key));
		idxmap[vertices[i]->key] = i;
	}
	for(int i=0; i<forrest.size(); ++i)
	{
		int k1 = idxmap[forrest[i]->u->key];
		int k2 = idxmap[forrest[i]->v->key];
		merge(nodes[k1], nodes[k2]);
	}
	vector<Node<Triangulation::_Internal::_edge*>*> cl = clusters(nodes);
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
	vector<Vertex<Triangulation::_Internal::_edge*>*> vertices;
	vector<Edge<Triangulation::_Internal::_edge*>*> edges;
	GraphFactory<Triangulation::_Internal::_edge*>* factory = GraphFactory<Triangulation::_Internal::_edge*>::GetInstance();
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(trmap.edges[i]));
	}
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* left = trebles[i].left;
		Triangulation::_Internal::_edge* right = trebles[i].right;
		int k1 = eidxmap[left];
		int k2 = eidxmap[right];
		Edge<Triangulation::_Internal::_edge*>* arc1 = factory->makeEdge(vertices[i], vertices[k1], 1.0-trebles[i].leftWeight);
		Edge<Triangulation::_Internal::_edge*>* arc2 = factory->makeEdge(vertices[i], vertices[k2], 1.0-trebles[i].rightWeight);
		vertices[i]->Add(arc1);
		edges.push_back(arc1);
		vertices[i]->Add(arc2);
		edges.push_back(arc2);
		Edge<Triangulation::_Internal::_edge*>* arc3 = factory->makeEdge(vertices[k1], vertices[i], 1.0-trebles[i].leftWeight);
		Edge<Triangulation::_Internal::_edge*>* arc4 = factory->makeEdge(vertices[k2], vertices[i], 1.0-trebles[i].rightWeight);
		vertices[k1]->Add(arc3);
		edges.push_back(arc3);
		vertices[k2]->Add(arc4);
		edges.push_back(arc4);
	}
	vector<Edge<Triangulation::_Internal::_edge*>*> mst = Kruskal(edges, vertices);
	vector<Edge<Triangulation::_Internal::_edge*>*> forrest;
	float maxW = 0;
	for(int i=0; i<mst.size(); ++i)
	{
		maxW = Max(maxW, mst[i]->w);
		if(mst[i]->w < thres)
		{
			forrest.push_back(mst[i]);
			//printf("%d: %d - %d -- %f\n", forrest.size(), pidxmap[mst[i]->u->key]+1, pidxmap[mst[i]->v->key]+1, mst[i]->w);
		}
	}
	printf("From %d edges, %d mst edges and %d forrest edges are kept with thres=%f\n", edges.size(), mst.size(), forrest.size(), thres);
	vector<int> label = connectedComponents(vertices, forrest);
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
		const int dims[] = {trmap.faces.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[1]->p));
			int j3 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[2]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.edges.size(), 4};
		vector<float> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], (float)(j1+1));
			SetData2(F, i, 1, dims[0], dims[1], (float)(j2+1));
			SetData2(F, i, 2, dims[0], dims[1], (float)label[i]);
			SetData2(F, i, 2, dims[0], dims[1], trebles[i].fitness);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		const int dims[] = {mst.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<mst.size(); ++i)
		{
			int j1 = eidxmap[mst[i]->u->key];
			int j2 = eidxmap[mst[i]->v->key];
			SetData2(F, i, 0, dims[0], dims[1], (j1+1));
			SetData2(F, i, 1, dims[0], dims[1], (j2+1));
			if(find(forrest.begin(), forrest.end(), mst[i])==forrest.end())
			{
				SetData2(F, i, 2, dims[0], dims[1], 0);
			}
			else
			{
				SetData2(F, i, 2, dims[0], dims[1], 1);
			}
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	factory->Clean();
	mexUnlock();
}

