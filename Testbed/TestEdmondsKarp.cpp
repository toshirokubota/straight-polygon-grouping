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
#include <Graph.h>
#include <GraphFactory.h>
#include <EdmondsKarp.h>
GraphFactory<string>* GraphFactory<string>::_instance = NULL;

namespace CSCI381
{
	/*
	Breadth First Search.
	*/
	AdjacencyMatrix
	BFS(int s, 
		AdjacencyMatrix& adj)
	{
		AdjacencyMatrix tree(adj.size());
		vector<float> d(adj.size(), numeric_limits<float>::infinity());
		vector<VertexColor> color(adj.size(), White);
		d[s] = 0;

		queue<int> q;
		q.push(s);
		while(q.empty() == false)
		{
			int n = q.front();
			q.pop();

			for(int i=0; i<adj.size(); ++i)
			{
				if(adj.Get(n, i)>0)
				{
					if(color[i] == White)
					{
						color[i] = Gray;
						d[i] = d[n] + 1;
						tree.Set(n, i, d[i]);
						q.push(i);
					}
				}
			}
			color[n] = Black;
		}
		return tree;
	}

	/*
	Given a network graph adn a flow, the function constructs a residual
	network and returns it as an AdjacentMatrix.
	*/
	AdjacencyMatrix
	GetResidual(AdjacencyMatrix& capacity,
		AdjacencyMatrix& flow)
	{
		AdjacencyMatrix residual = capacity;
		for(int i=0; i<flow.size(); ++i)
		{
			for(int j=0; j<flow.size(); ++j)
			{
				if(flow.Get(i, j) > 0)
				{
					residual.Set(i, j, capacity.Get(i, j) - flow.Get(i, j));
					residual.Set(j, i, capacity.Get(j, i) + flow.Get(i, j));
				}
			}
		}
		return residual;
	}

	/*
	This function finds an augmenting path from a residual network.  It 
	returns true if an augmenting path is found, and returns false otherwise.
	If found, the path is stored in an AdjacentMatrix, path.
	*/
	bool
	AugmentingPath(AdjacencyMatrix& path,
		AdjacencyMatrix& residual,
		int s,
		int t)
	{
		AdjacencyMatrix tree = BFS(s, residual);
		while(t != s)
		{
			int id = -1;
			for(int i=0; i<tree.size(); ++i)
			{
				if(tree.Get(i, t) > 0)
				{
					id = i;
					break;
				}
			}
			if(id < 0)
			{
				break;
			}
			else
			{
				path.Set(id, t, residual.Get(id, t));
				t = id;
			}
		}
		if(t != s)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
#include <cassert>

	/*
	The function adds an additional flow along the given path to the current flow.
	*/
	void
	updateFlow(AdjacencyMatrix& flow,
		AdjacencyMatrix& path,
		int s, 
		int t)
	{
		//find the minimum capacity along the augumenting path
		float inc = numeric_limits<float>::infinity();
		for(int i=0; i<path.size(); ++i)
		{
			for(int j=0; j<path.size(); ++j)
			{
				if(path.Get(i, j) > 0)
				{
					inc = min(inc, path.Get(i,j));
				}
			}
		}

		//add the minimum capacity to the current flow along the augmenting path
		for(int i=0; i<path.size(); ++i)
		{
			for(int j=0; j<path.size(); ++j)
			{
				if(path.Get(i, j) > 0)
				{
					flow.Set(i, j, flow.Get(i, j) + inc);
					//flow.Set(j, i, -flow.Get(i, j));
				}
			}
		}
	}

	AdjacencyMatrix
	EdmondsKarp(AdjacencyMatrix& adj,
		int s,
		int t)
	{
		AdjacencyMatrix flow(adj.size());
		adj.Print();

		do
		{
			AdjacencyMatrix residual = GetResidual(adj, flow);
			cout << "Flow:" << endl; flow.Print();
			cout << "Residual:" << endl; residual.Print();
			AdjacencyMatrix path(residual.size());
			if(AugmentingPath(path, residual, s, t))
			{
				cout << "Path:" << endl; path.Print();
				updateFlow(flow, path, s, t);
			}
			else
			{
				break;
			}
		}
		while(1);

		return flow;
	}
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: T = Testbed(infile)");
		return;
	}

	ifstream in("C:\\Toshiro\\Classes\\2012\\CSCI381\\Programs\\EdmondsKarp\\data3.txt");
	if(in.fail() == false)
	{
		vector<Vertex<string>*> vnodes;
		AdjacencyMatrix adj;
		if(ReadFromFile(in, vnodes, adj))
		{
			for(int i=0; i<vnodes.size(); ++i)
			{
				printf("%s: ", vnodes[i]->key.c_str());
				for(int j=0; j<vnodes[i]->aList.size(); ++j)
				{
					Edge<string>* e = vnodes[i]->aList[j];
					printf("(%s,%s,%f), ", e->u->key.c_str(), e->v->key.c_str(), e->w);
				}
				printf("\n");
			}
			//EdmondsKarp(adj, 0, vnodes.size()-1);
			float total = EdmondsKarp(vnodes, vnodes[0], vnodes[vnodes.size()-1]);
			vector<bool> label = minCut(vnodes, vnodes[0]);
			printf("Total flow = %f\n", total);
			printf("Cut:\n");
			for(int i=0; i<vnodes.size(); ++i)
			{
				if(label[i]) printf("\t%s\n", vnodes[i]->key.c_str());
				else printf("%s\n", vnodes[i]->key.c_str());
			}
		}
	}
	mexUnlock();
}

