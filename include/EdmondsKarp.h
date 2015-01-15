// EdmondsKarp.cpp : Defines the entry point for the console application.
//

#include <Graph.h>
#include <queue>
using namespace std;
/*
Breadth First Search.
*/
template <class T>
float
BFS(Vertex<T>* node, 
	vector<Vertex<T>*>& vnodes)
{
	//set the color, d, and pi to White, infinity, and NULL
	for(int i=0; i<vnodes.size(); ++i)
	{
		vnodes[i]->Reset();
	}

	node->d = 0;
	//enqueue the pointer to the node
	queue<Vertex<T>*> q;
	q.push(node);
	float min_weight = std::numeric_limits<float>::infinity();
	while(q.empty() == false)
	{
		//This is what you need to implement
		Vertex<T>* u = q.front(); 
		q.pop();
		for(int i=0; i<u->aList.size(); ++i)
		{
			if(u->aList[i]->w > 0) //use only edges with non-zero weights
			{
				Vertex<T>* v = u->aList[i]->v;
				if(v->color == White)
				{
					v->color = Gray;
					v->d = u->d + 1;
					v->pi = u;
					q.push(v);
				}
				min_weight = Min(min_weight, u->aList[i]->w);
			}
		}
		u->color = Black;
	}
	return min_weight;
}

template<class T>
float
getResidualFlow(Vertex<T>* source, Vertex<T>* target)
{
	float min_weight = std::numeric_limits<float>::infinity();
	Vertex<T>* q = target;
	while(q->pi)
	{
		Vertex<T>* u = q->pi;
		Edge<T>* e = NULL;
		for(int i=0; i<u->aList.size(); ++i)
		{
			if(u->aList[i]->v == q)
			{
				e = u->aList[i];
				break;
			}
		}
		assert(e);
		if(e->w < min_weight) min_weight = e->w;

		q = q->pi;
	}
	assert(q == source);
	return min_weight;
}

template <class T>
void
updateFlowNetwork(Vertex<T>* target, float residual_flow)
{
	Vertex<T>* q = target;
	while(q->pi)
	{
		{ //forward flow
			Vertex<T>* u = q->pi;
			FlowEdge<T>* e = NULL;
			for(int i=0; i<u->aList.size(); ++i)
			{
				if(u->aList[i]->v == q)
				{
					e = (FlowEdge<T>*) u->aList[i];
					break;
				}
			}
			assert(e);
			e->w -= residual_flow;
			e->flow += residual_flow;
		}
		{ //reverse flow
			Vertex<T>* u = q;
			FlowEdge<T>* e = NULL;
			for(int i=0; i<u->aList.size(); ++i)
			{
				if(u->aList[i]->v == q->pi)
				{
					e = (FlowEdge<T>*) u->aList[i];
					break;
				}
			}
			//assert(e);
			if(e)
			{
				e->w += residual_flow;
			}
		}
		q = q->pi;
	}
}

template <class T>
float
EdmondsKarp(vector<Vertex<T>*> vnodes, Vertex<T>* source, Vertex<T>* target)
{
	float total_flow = 0;
	int iter = 0;
	do
	{
		iter++;
		float min_weight = BFS(source, vnodes);
		if(target->pi == NULL) break;
		float residual_flow = getResidualFlow(source, target);
		total_flow += residual_flow;
		updateFlowNetwork(target, residual_flow);
	}
	while(1);

	return total_flow;
}

#include <map>
#include <set>
using namespace std;
#include <mex.h>

template <class T>
vector<bool> 
minCut(vector<Vertex<T>*> vnodes, Vertex<T>* source)
{
	vector<bool> cluster(vnodes.size(), false);
	map<Vertex<T>*,int> indexmap;
	for(int i=0; i<vnodes.size(); ++i)
	{
		indexmap[vnodes[i]] = i;
	}
	set<Vertex<T>*> Q;
	Q.insert(source);
	while(Q.empty() == false)
	{
		set<Vertex<T>*> Q0 = Q;
		Q.clear();
		for(set<Vertex<T>*>::iterator it=Q0.begin(); it!=Q0.end(); it++)
		{
			Vertex<T>* u = *it;
			int k = indexmap[u];
			cluster[k] = true;
			for(int j=0; j<u->aList.size(); ++j)
			{
				Edge<T>* e = u->aList[j];
				if(e->w > 0)
				{
					int k2 = indexmap[e->v];
					if(cluster[k2]==false)
					{
						Q.insert(e->v);
					}
				}
			}
		}
	}

	return cluster;
}

