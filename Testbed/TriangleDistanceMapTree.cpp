#include <iostream>
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
#include <szDistanceTransform.h>
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>
#include <intersectionConvexPolygons.h>
#include <szConvexHull2D.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Tree.h>

typedef Triangulation::_Internal::_triangle* GraphNodeType; 
GraphFactory<GraphNodeType>* GraphFactory<GraphNodeType>::_instance = 0;

void
trace(Vertex<GraphNodeType>* node, float weight)
{
	for(int i=0; i<node->aList.size(); ++i)
	{
		Vertex<GraphNodeType>* next = node->aList[i]->v;
		if(next->pi == node)
		{
			trace(next, weight);
		}
	}
	float sum = 0;
	int count = 0;
	float minval = std::numeric_limits<float>::infinity();
	for(int i=0; i<node->aList.size(); ++i)
	{
		Vertex<GraphNodeType>* next = node->aList[i]->v;
		if(next->pi == node)
		{
			sum += next->f;
			count ++;
			minval = Min(minval, next->f);
		}
	}
	if(count > 0)
	{
		node->f = (node->f + weight * sum / count) / (1.0f + weight);
		if(node->f > minval)
		{
			node->f = minval; //make sure the value is not larger than its children.
		}
	}
}

/*
Trace a tree in a breadth-first, then take a weighted average of the distance between itself and its children.
*/
void
smoothenDistance(vector<Vertex<GraphNodeType>*>& nodes, Vertex<GraphNodeType>* root, float weight)
{
	for(int i=0; i<nodes.size(); ++i)
	{
		nodes[i]->f = nodes[i]->d;
	}
	trace(root, weight);
	float offset = root->f;
	for(int i=0; i<nodes.size(); ++i)
	{
		nodes[i]->f -= offset;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 4 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [J] = Testbed(points, faces, fitness, seed)");
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
	//triangle (indices to the points)
	vector<Triangulation::_Internal::_indexed_triangle> T;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
		T = indices2structs(points, T0, dimsT);
	}
	//boundary
	vector<float> boundary;
	{
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(boundary, prhs[2], classIdT, ndimT, &dimsT);
	}
	int seed = -1; //seed triangle
	{
		mxClassID classMode;
		ReadScalar(seed,prhs[3],classMode);
		seed--; //change from MATLAB index.
	} 
	float thres = 10.0;
	if(nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[4],classMode);
	} 
	float weight = 2.0;
	if(nrhs >= 6)
	{
		mxClassID classMode;
		ReadScalar(weight,prhs[5],classMode);
	} 

	Triangulation::Triangulator trmap(points, T);
	Triangulation::_Internal::_triangle* start = trmap.faces[seed];

	map<Triangulation::_Internal::_vertex*,int> idxmap;
	for(int i=0; i<trmap.points.size(); ++i)
	{
		idxmap[trmap.points[i]] = i;
	}
	map<Triangulation::_Internal::_triangle*,int> fmap;
	map<Triangulation::_Internal::_triangle*,CParticleF> center;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		Triangulation::_Internal::_triangle* tr = trmap.faces[i];
		fmap[tr] = i;
		center[trmap.faces[i]] = CParticleF((tr->vertices[0]->p.m_X+tr->vertices[1]->p.m_X+tr->vertices[2]->p.m_X)/3.0,
			(tr->vertices[0]->p.m_Y+tr->vertices[1]->p.m_Y+tr->vertices[2]->p.m_Y)/3.0);
	}
	map<Triangulation::_Internal::_edge*,int> edxmap;
	map<Triangulation::_Internal::_edge*,float> bdrmap;
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		edxmap[trmap.edges[i]] = i;
		bdrmap[trmap.edges[i]] = boundary[i];
	}

	vector<Vertex<GraphNodeType>*> nodes(trmap.faces.size());
	GraphFactory<GraphNodeType>* factory = GraphFactory<GraphNodeType>::GetInstance();
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		nodes[i] = factory->makeVertex(trmap.faces[i]);
		nodes[i]->d = std::numeric_limits<float>::infinity();
	}
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		Triangulation::_Internal::_triangle* tr = trmap.faces[i];
		for(int j=0; j<3; ++j)
		{
			Triangulation::_Internal::_edge* ed = tr->edges[j];
			Triangulation::_Internal::_triangle* tr2 = ed->faces[0]==tr ? ed->faces[1]: ed->faces[0];
			if(tr2)
			{
				int k = fmap[tr2];
				nodes[i]->Add(factory->makeEdge(nodes[i], nodes[k], 1.0));
			}
		}
	}
	nodes[fmap[start]]->d = 0;
	int tic=0;
	while(true)
	{
		tic++;
		bool bChanged = false;
		for(int i=0; i<nodes.size(); ++i)
		{
			Triangulation::_Internal::_triangle* tr = nodes[i]->key;
			CParticleF c1 = center[tr];
			float d = nodes[i]->d;
			for(int j=0; j<nodes[i]->aList.size(); ++j)
			{
				Vertex<GraphNodeType>* n2 = nodes[i]->aList[j]->v;
				float d2 = n2->d;
				Triangulation::_Internal::_triangle* tr2 = n2->key;
				CParticleF c2 = center[tr2];
				float sep = Distance(c1, c2);
				Triangulation::_Internal::_edge* ed = NULL;
				for(int k=0; k<3; ++k)
				{
					if((tr->edges[k]->faces[0]==tr && tr->edges[k]->faces[1]==tr2) || (tr->edges[k]->faces[1]==tr && tr->edges[k]->faces[0]==tr2))
					{
						ed = tr->edges[k];
						break;
					}
				}
				if(bdrmap[ed] <= 0)
				{
					float d3 = d + sep;
					if(d2 > d3)
					{
						n2->d = d3;
						n2->pi = nodes[i];
						bChanged = true;
					}
				}
			}
		}
		if(bChanged == false) break;
	}

	smoothenDistance(nodes, nodes[fmap[start]], weight);
	
	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], trmap.points[i]->p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], trmap.points[i]->p.m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trmap.edges.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = idxmap[trmap.edges[i]->vertices[0]];
			int j2 = idxmap[trmap.edges[i]->vertices[1]];
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.faces.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = idxmap[trmap.faces[i]->vertices[0]];
			int j2 = idxmap[trmap.faces[i]->vertices[1]];
			int j3 = idxmap[trmap.faces[i]->vertices[2]];
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		const int dims[] = {trmap.faces.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			float d = nodes[i]->d;
			float f = nodes[i]->f;
			if(d == std::numeric_limits<float>::infinity())
			{
				SetData2(F, i, 0, dims[0], dims[1], -1.0f);
			}
			else
			{
				SetData2(F, i, 0, dims[0], dims[1], d);
			}
			if(f == std::numeric_limits<float>::infinity())
			{
				SetData2(F, i, 1, dims[0], dims[1], -1.0f);
			}
			else
			{
				SetData2(F, i, 1, dims[0], dims[1], f);
			}
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

