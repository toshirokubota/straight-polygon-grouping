#ifndef ___GRAPH_CSCI381_H___
#define ___GRAPH_CSCI381_H___
#include <vector>
#include <limits>
#include <iostream>

enum VertexColor {Black, White, Gray, Red, Green, Blue};
enum EdgeType{Tr, Back, Forward, Cross};

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

	void Add(Edge<T>* edge)
	{
		aList.push_back(edge);
	}

	Edge<T>* findEdge(Vertex<T>* v)
	{
		for(int i=0; i<aList.size(); ++i)
		{
			if(aList[i]->v == v)
			{
				return aList[i];
			}
		}
		return NULL;
	}

	void Reset()
	{
		color = White;
		d = numeric_limits<float>::infinity();
		f = numeric_limits<float>::infinity();
		pi = NULL;
		//children.clear();
	}
	bool operator ==(Vertex<T>* n)
	{
		return key == n->key;
	}

	T key;
	std::vector<Edge<T>*> aList;
	//std::vector<Vertex<T>*> children;
	VertexColor color;
	float d;
	float f;
	Vertex<T>* pi;
};

template <class T>
class Edge
{
public:
	Edge(Vertex<T>* u=0, Vertex<T>* v=0, float w=0, EdgeType type = Tr)
	{
		this->u = u;
		this->v = v;
		this->w = w;
		this->type = type;
	}
	bool operator ==(Edge<T>* n)
	{
		return u == n->u && v == n->v;
	}
	bool operator <(Edge<T>* n)
	{
		return w < n->w;
	}
	virtual ~Edge() {}
	Vertex<T>* u;
	Vertex<T>* v;
	float w;
	EdgeType type;
};

template <class T>
class FlowEdge: public Edge<T>
{
public:
	FlowEdge(Vertex<T>* u=0, Vertex<T>* v=0, float w=0, EdgeType type = Tr): Edge<T>(u, v, w, type)
	{
		this->flow = 0;
	}
	virtual ~FlowEdge() {}
	float flow;
};

class AdjacencyMatrix
{
public:
	AdjacencyMatrix(int n=0): weights(n*n, 0)
	{
		N = n;
	}
	void Set(int i, int j, float w = 1)
	{
		weights[i*N+j] = w;
	}
	void SetSymmetric(int i, int j, float w = 1)
	{
		Set(i, j, w);
		Set(j, i, w);
	}
	float Get(int i, int j) const 
	{
		return weights[i*N+j];
	}
	int size() const {return N;}
	void Print() const
	{
		printf("Adjacency Matrix:\n");
		for(int i=0; i<N; ++i)
		{
			for(int j=0; j<N; ++j)
			{
				printf("%f ", Get(i,j));
			}
			printf("\n");
		}
	}
	std::vector<float> GetWeights() {return weights;}

protected:
	int N;
	std::vector<float> weights;
};

class AdjacencyPredecessorMatrix: public AdjacencyMatrix
{
public:
	AdjacencyPredecessorMatrix(int n=0): AdjacencyMatrix(n), pi(n*n, -1)
	{
	}
	AdjacencyPredecessorMatrix(const AdjacencyMatrix& M)
	{
		N = M.size();
		weights = std::vector<float>(N*N);
		for(int i=0; i<N; ++i)
		{
			for(int j=0; j<N; ++j)
			{
				Set(i, j, M.Get(i, j));
			}
		}
		pi = std::vector<int>(N*N, -1);
	}
	AdjacencyPredecessorMatrix& operator =(const AdjacencyMatrix& M)
	{
		N = M.size();
		weights = std::vector<float>(N*N);
		for(int i=0; i<N; ++i)
		{
			for(int j=0; j<N; ++j)
			{
				Set(i, j, M.Get(i, j));
			}
		}
		pi = std::vector<int>(N*N, -1);
		return *this;
	}
	void SetPi(int i, int j, int p)
	{
		pi[i*N+j] = p;
	}
	int GetPi(int i, int j) const 
	{
		return pi[i*N+j];
	}

protected:
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

	int nE;
	in >> nE;
	adj = AdjacencyMatrix(N);
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
		Edge<T>* edge = new FlowEdge<T>(vnodes[m], vnodes[n], w);
		vnodes[m]->Add(edge);
		vedges.push_back(edge);

		adj.Set(m, n, w);
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
std::vector<Vertex<T>*>
DFS(int k,
	const AdjacencyMatrix& adj,
	int& time,
	std::vector<Vertex<T>*> vnodes)
{
	time++;
	vnodes[k]->color = Gray;
	vnodes[k]->d = time;
	std::vector<Vertex<T>*> trace;
	for(int i=0; i<adj.size(); ++i)
	{
		if(adj.Get(k, i)>0)
		{
			Vertex<T>* node = vnodes[i];
			if(node->color == White)
			{
				node->pi = vnodes[k];
				std::vector<Vertex<T>*> trace0 = DFS(i, adj, time, vnodes);
				trace.insert(trace.end(), trace0.begin(), trace0.end());
			}
		}
	}
	vnodes[k]->color = Black;
	vnodes[k]->f = ++time;
	trace.push_back(vnodes[k]);
	return trace;
}

template<class T>
int
partitionEdges(std::vector<Edge<T>*>& edges, int a, int b)
{
	Edge<T>* pivot = edges[b];
	int i=a;
	for(int j=a; j<b; ++j)
	{
		if(edges[j]->w < pivot->w)
		{
			swap(edges[i], edges[j]);
			i++;
		}
	}
	swap(edges[b], edges[i]);
	return i;
}

template<class T>
void
sortEdges(std::vector<Edge<T>*>& edges, int i, int j)
{
	if(i>=j)
	{
		return;
	}
	int k = partitionEdges(edges, i, j);
	sortEdges(edges, i, k-1);
	sortEdges(edges, k+1, j);
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

template <class T>
Vertex<T>*
findVertex(const std::vector<Vertex<T>*>& vertices,
		   const T& t)
{
	for(int i=0; i<vertices.size(); ++i)
	{
		if(vertices[i]->key == t)
		{
			return vertices[i];
		}
	}
	return NULL;
}

#endif /* ___GRAPH_CSCI381_H___*/