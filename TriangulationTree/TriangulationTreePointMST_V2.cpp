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
#include <cassert>

GraphFactory<Triangulation::_Internal::_vertex*>* GraphFactory<Triangulation::_Internal::_vertex*>::_instance = NULL;

vector<int> 
connectedComponents(vector<Vertex<Triangulation::_Internal::_vertex*>*>& vertices, 
					vector<Edge<Triangulation::_Internal::_vertex*>*>& forrest) 
{
	map<Triangulation::_Internal::_vertex*,int> idxmap;
	vector<Node<Triangulation::_Internal::_vertex*>*> nodes;
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
	vector<Node<Triangulation::_Internal::_vertex*>*> cl = clusters(nodes);
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

	Triangulation::Triangulator trmap(points);
	vector<float> W(trmap.edges.size());
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		W[i] = trmap.edges[i]->Length();
	}

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
	/*map<Triangulation::_Internal::_triangle*,int> fidxmap;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		fidxmap[trmap.faces[i]]=i;
	}*/

	vector<Vertex<Triangulation::_Internal::_vertex*>*> vertices;
	vector<Edge<Triangulation::_Internal::_vertex*>*> edges;
	GraphFactory<Triangulation::_Internal::_vertex*>* factory = GraphFactory<Triangulation::_Internal::_vertex*>::GetInstance();
	for(int i=0; i<trmap.points.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(trmap.points[i]));
	}
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_vertex* pu = trmap.edges[i]->vertices[0];
		Triangulation::_Internal::_vertex* pv = trmap.edges[i]->vertices[1];
		int k1 = pidxmap[pu];
		int k2 = pidxmap[pv];
		Edge<Triangulation::_Internal::_vertex*>* arc1 = factory->makeEdge(vertices[k1], vertices[k2], W[i]);
		Edge<Triangulation::_Internal::_vertex*>* arc2 = factory->makeEdge(vertices[k2], vertices[k1], W[i]);
		vertices[k1]->Add(arc1);
		edges.push_back(arc1);
		vertices[k2]->Add(arc2);
		edges.push_back(arc2);
	}
	vector<Edge<Triangulation::_Internal::_vertex*>*> mst = Kruskal(edges, vertices);
	vector<bool> vbflag(trmap.edges.size(), false);
	for(int i=0; i<mst.size(); ++i)
	{
		Triangulation::_Internal::_vertex* u = mst[i]->u->key;
		Triangulation::_Internal::_vertex* v = mst[i]->v->key;
		Triangulation::_Internal::_edge* ed = findEdge(u, v);
		vbflag[eidxmap[ed]] = true;
	}

	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 3};
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
			int j1 = pidxmap[trmap.faces[i]->vertices[0]];
			int j2 = pidxmap[trmap.faces[i]->vertices[1]];
			int j3 = pidxmap[trmap.faces[i]->vertices[2]];
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.edges.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], (j1+1));
			SetData2(F, i, 1, dims[0], dims[1], (j2+1));
			SetData2(F, i, 2, dims[0], dims[1], vbflag[i] ? 1: 0);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	factory->Clean();
	mexUnlock();
}

