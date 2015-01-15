#ifndef ___GRAPH_CSCI381_H___
#define ___GRAPH_CSCI381_H___
#include <vector>
#include <limits>
#include <iostream>

enum VertexColor {Black, White, Gray};
enum EdgeType{Tree, Back, Forward, Cross};

template<class T>
class Edge;

template <class T>
class Vertex
{
public:
	Vertex(T key)
	{
		this->key = key;
		color = White;
		d = 0;
		f = 0;
		pi = NULL;
	}
	~Vertex()
	{
		for(int i=0; i<aList.size(); ++i)
		{
			delete aList[i];
		}
	}

	void Add(Edge<T>* edge)
	{
		aList.push_back(edge);
	}

	void Reset()
	{
		color = White;
		d = numeric_limits<float>::infinity();
		pi = NULL;
	}

	T key;
	std::vector<Edge<T>*> aList;
	VertexColor color;
	float d;
	float f;
	Vertex<T>* pi;
};

template <class T>
class Edge
{
public:
	Edge(Vertex<T>* u=0, Vertex<T>* v=0, float w=0, float f=0, EdgeType type = Tree)
	{
		this->u = u;
		this->v = v;
		this->w = w;
		this->type = type;
		this->flow = f;
	}
	Vertex<T>* u;
	Vertex<T>* v;
	float w;
	float flow;
	EdgeType type;
};

class AdjacencyMatrix
{
public:
	AdjacencyMatrix(int n=0): weights(n*n, std::numeric_limits<float>::infinity()), pi(n*n, -1)
	{
		N = n;
	}
	void Set(int i, int j, float w = 1)
	{
		weights[i*N+j] = w;
	}
	void SetPi(int i, int j, int p)
	{
		pi[i*N+j] = p;
	}
	float Get(int i, int j) const 
	{
		return weights[i*N+j];
	}
	int GetPi(int i, int j) const 
	{
		return pi[i*N+j];
	}
	int size() const {return N;}
	void Print() const
	{
		std::cout << "Adjacency Weight Matrix:" << std::endl;
		for(int i=0; i<N; ++i)
		{
			for(int j=0; j<N; ++j)
			{
				std::cout << Get(i,j) << "\t";
			}
			std::cout << std::endl;
		}
	}
	void PrintPi() const
	{
		std::cout << "Adjacency PI Matrix:" << std::endl;
		for(int i=0; i<N; ++i)
		{
			for(int j=0; j<N; ++j)
			{
				std::cout << GetPi(i,j) << "\t";
			}
			std::cout << std::endl;
		}
	}
protected:
	int N;
	std::vector<float> weights;
	std::vector<int> pi;
};

#include <fstream>
#include <algorithm>
template<class T>
bool
ReadFromFile(std::istream& in,
			 std::vector<Vertex<T>*>& vnodes,
			 AdjacencyMatrix& adj)
{
	int N;
	in >> N;
	vnodes = std::vector<Vertex<T>*>(N);
	std::vector<T> vkeys(N);
	for(int i=0; i<N; ++i)
	{
		T key;
		in >> key;
		vnodes[i] = new Vertex<T>(key);
		vkeys[i] = key;
	}
	if(in.bad())
	{
		std::cerr << "Error while reading a vertex list." << std::endl;
		std::cerr << "The list of vertices are: " << std::endl;
		for(int i=0; i<vnodes.size(); ++i)
		{
			std::cerr << vnodes[i]->key << std::endl;
		}
		return false;
	}
	adj = AdjacencyMatrix(N);
	for(int i=0; i<N; ++i) //diagonal of the adjacency matrix to 0. others are set to infinity
	{
		adj.Set(i, i, 0);
		adj.SetPi(i, i, i);
	}

	int nE;
	in >> nE;
	std::vector<Edge<T>*> vedges;
	bool bSuccess = true;
	int count = 0;

	for(int i=0; i<nE; ++i, ++count)
	{
		T s, t;
		float w;
		in >> s >> t >> w;

		std::vector<T>::iterator ps = std::find(vkeys.begin(), vkeys.end(), s);
		std::vector<T>::iterator pt = std::find(vkeys.begin(), vkeys.end(), t);
		if(ps == vkeys.end())
		{
			std::cerr << "Could not find " << s << " in the vertex list." << std::endl;
			bSuccess = false;
			break;
		}
		if(pt==vkeys.end())
		{
			std::cerr << "Could not find " << t << " in the vertex list." << std::endl;
			bSuccess = false;
			break;
		}
		int m = (int)(ps - vkeys.begin());
		int n = (int)(pt - vkeys.begin());
		Edge<T>* edge = new Edge<T>(vnodes[m], vnodes[n], w);
		vnodes[m]->Add(edge);
		vedges.push_back(edge);

		adj.Set(m, n, w);
		adj.SetPi(m, n, m);
	}
	if(!bSuccess || in.bad())
	{
		std::cerr << "Error while reading an edge list." << std::endl;
		std::cerr << "# of edges read is " << count << std::endl;
		std::cerr << "The list of edges are: " << std::endl;
		for(int i=0; i<vedges.size(); ++i)
		{
			Edge<T>* edge = vedges[i];
			std::cerr << edge->u->key << " " << edge->v->key << " " << edge->w << " " << std::endl;
		}
		return false;
	}
	else
	{
		return true;
	}
}

template <class T>
std::vector<Edge<T>*>
GetAllEdges(const std::vector<Vertex<T>*>& vnodes)
{
	std::vector<Edge<T>*> edges;
	for(int i=0; i<vnodes.size(); ++i)
	{
		for(int j=0; j<vnodes[i]->aList.size(); ++j)
		{
			edges.push_back(vnodes[i]->aList[j]);
		}
	}
	return edges;
}

template <class T>
std::vector<Vertex<T>*>
GetAllVertices(const std::vector<Edge<T>*>& edges)
{
	std::vector<Vertex<T>*> nodes;
	for(int i=0; i<edges.size(); ++i)
	{
		if(std::find(nodes.begin(), nodes.end(), edges[i]->u) == nodes.end())
		{
			nodes.push_back(edges[i]->u);
		}
		if(std::find(nodes.begin(), nodes.end(), edges[i]->v) == nodes.end())
		{
			nodes.push_back(edges[i]->v);
		}
	}
	return nodes;
}


#endif /* ___GRAPH_CSCI381_H___*/