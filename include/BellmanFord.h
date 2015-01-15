// BellmanFord.cpp : Defines the entry point for the console application.
//

#include <Graph.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

template <class T>
bool
negativeCycle(const vector<Edge<T>*>& edges)
{
	for(int i=0; i<edges.size(); ++i)
	{
		Edge<T>* e = edges[i];
		if(e->v->d > e->u->d + e->w)
		{
			return false;
		}
	}
	return true;
}

template <class T>
bool
BellmanFord(vector<Vertex<T>*>& vertices,
			vector<Edge<T>*>& edges,
			Vertex<T>* source)
{
	for(int i=0; i<vertices.size(); ++i)
	{
		vertices[i]->Reset();
	}
	source->d = 0;
	source->pi = NULL;

	for(int i=0; i<vertices.size() - 1; ++i)
	{
		for(int j=0; j<edges.size(); ++j)
		{
			Edge<T>* e = edges[j];
			if(e->v->d > e->u->d + e->w)
			{
				e->v->d = e->u->d + e->w;
				e->v->pi = e->u;
			}
		}
	}

	return negativeCycle(edges);
}
