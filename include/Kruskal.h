// Kruskal.cpp : Defines the entry point for the console application.
//

#include <vector>
#include <algorithm>
using namespace std;
#include <Graph.h>
#include <DisjointSet.h>

template<class T>
vector<Edge<T>*>
Kruskal(vector<Edge<T>*>& edges, 
		vector<Vertex<T>*>& vertices)
{
	vector<Edge<T>*> mst;
	sortEdges(edges, 0, edges.size()-1);

	vector<Node<Vertex<T>*>*> nodes;
	for(int i=0; i<vertices.size(); ++i)
	{
		nodes.push_back(makeset(vertices[i]));
	}

	for(int i=0; i<edges.size(); ++i)
	{
		int ui = distance(vertices.begin(), find(vertices.begin(), vertices.end(), edges[i]->u));
		int uj = distance(vertices.begin(), find(vertices.begin(), vertices.end(), edges[i]->v));
		if(findset(nodes[ui]) != findset(nodes[uj]))
		{
			merge(nodes[ui], nodes[uj]);
			mst.push_back(edges[i]);
		}
	}
	for(int i=0; i<vertices.size(); ++i)
	{
		delete nodes[i];
	}
	return mst;
}

#include <set>
template<class T>
vector<Edge<T>*>
Kruskal(vector<Edge<T>*>& edges)
{
	set<Vertex<T>*> V;
	for(int i=0; i<edges.size(); ++i)
	{
		V.insert(edges[i]->u);
		V.insert(edges[i]->v);
	}
	vector<Vertex<T>*> vertices;
	for(set<Vertex<T>*>::iterator p = V.begin(); p!=V.end(); p++)
	{
		vertices.push_back(*p);
	}
	return Kruskal(edges, vertices);
}

