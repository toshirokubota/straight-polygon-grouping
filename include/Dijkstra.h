#ifndef ___DIJKSTRA_H___
#define ___DIJKSTRA_H___
#include <Graph.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

/*
Brute force O(n) implementation of extracting a vertex with the minimum d.
*/
template <class T>
Vertex<T>*
extractMinimum(std::vector<Vertex<T>*>& Q)
{
	float minD = Q[0]->d;
	int minid = 0;
	for(int i=1; i<Q.size(); ++i)
	{
		if(minD > Q[i]->d)
		{
			minD = Q[i]->d;
			minid = i;
		}
	}
	Vertex<T>* minV = Q[minid];
	Q.erase(Q.begin()+minid);
	return minV;
}

template<class T>
void
Dijkstra(std::vector<Vertex<T>*>& vertices,
		 Vertex<T>* source)
{
	for(int i=0; i<vertices.size(); ++i)
	{
		vertices[i]->Reset();
	}
	source->d = 0;

	vector<Vertex<T>*> Q = vertices;
	while(Q.empty() == false)
	{
		Vertex<T>* u = extractMinimum(Q);
		for(int i=0; i<u->aList.size(); ++i)
		{
			Edge<T>* e = u->aList[i];
			Vertex<T>* v = e->v;
			if(v->d > u->d + e->w)
			{
				v->d = u->d + e->w;
				v->pi = u;
			}
		}
	}
}

/*
Dijkstra with a set of edges to choose from (instead of using aList).
*/
template<class T>
void
Dijkstra(std::vector<Vertex<T>*>& vertices,
		 std::vector<Edge<T>*>& edges,
		 Vertex<T>* source)
{
	for(int i=0; i<vertices.size(); ++i)
	{
		vertices[i]->Reset();
	}
	source->d = 0;

	vector<Vertex<T>*> Q = vertices;
	while(Q.empty() == false)
	{
		Vertex<T>* u = extractMinimum(Q);
		for(int i=0; i<edges.size(); ++i)
		{
			if(edges[i]->u == u)
			{
				Edge<T>* e = edges[i];
				Vertex<T>* v = e->v;
				if(v->d > u->d + e->w)
				{
					v->d = u->d + e->w;
					v->pi = u;
				}
			}
		}
	}
}

template <class T>
std::vector<Vertex<T>*>
tracePath(Vertex<T>* v, Vertex<T>* u=NULL)
{
	
	std::vector<Vertex<T>*> path;
	while(v && v!=u)
	{
		/*if(find(path.begin(), path.end(), v) != path.end())
		{
			v+=0;
		}*/
		path.insert(path.begin(), v);
		v = v->pi;
	}
	if(v) path.insert(path.begin(), v);
	return path;
}

#endif /* ___DIJKSTRA_H___ */