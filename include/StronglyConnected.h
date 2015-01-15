#ifndef ___STRONGLY_CONNECTED_H___
#define ___STRONGLY_CONNECTED_H___

#include <vector>
using namespace std;
#include <Graph.h>

template <class T>
vector<int>
StronglyConnectedComponents(vector<Vertex<T>*>& vnodes,
							const AdjacencyMatrix& adj)
{
	//Step 1: call DFS to compute finishing times for each vertex
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i]->color = White;
	}
	int time = 0;
	for(int i=0; i<vnodes.size(); ++i)
	{
		if(vnodes[i]->color == White)
		{
			DFS(i, adj, time, vnodes);
		}
	}
	//Step 2: compute G^T by transposing the adjacency matrix
	AdjacencyMatrix adjT(adj.size());
	for(int i=0; i<adj.size(); ++i)
	{
		for(int j=0; j<adj.size(); ++j)
		{
			adjT.Set(i, j, adj.Get(j, i));
		}
	}
	//Step 3: call DFS on G^T.  Consider the vertices in order of decreasing finishing time
	vector<pair<int,int>> pairs;
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i]->color = White;
		pairs.push_back(pair<int,int>(vnodes[i]->f, i));
	}
	sort(pairs.begin(), pairs.end());

	time = 0;
	int count = 1;
	
	vector<int> components(vnodes.size(), -1);
	int label = 1;
	for(int i=0; i<vnodes.size(); ++i)
	{
		//find a white vertex with the largest finishing time
		int index = pairs[vnodes.size()-1-i].second;
		if(vnodes[index]->color == White)
		{
			vector<Vertex<T>*> trace = DFS(index, adjT, time, vnodes);
			for(int k=0; k<trace.size(); ++k)
			{
				int l = distance(vnodes.begin(), find(vnodes.begin(), vnodes.end(), trace[k]));
				if(l >= 0 && l<vnodes.size())
				{
					components[l] = label;
				}
			}
			label ++;
		}
	}
	return components;
}

template <class T>
vector<int>
StronglyConnectedComponents(vector<Vertex<T>*>& vnodes)
{
	AdjacencyMatrix adj(vnodes.size());
	for(int i=0; i<vnodes.size(); ++i)
	{
		for(int j=0; j<vnodes[i]->aList.size(); ++j)
		{
			int k = distance(vnodes.begin(), find(vnodes.begin(), vnodes.end(), vnodes[i]->aList[j]->v));
			adj.Set(i, k, vnodes[i]->aList[j]->w);
		}
	}
	return StronglyConnectedComponents(vnodes, adj);
}

#endif /*  ___STRONGLY_CONNECTED_H___ */