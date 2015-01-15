#ifndef ___EDGE_PAIR_H___
#define ___EDGE_PAIR_H___
#include <Triangulation.h>

struct EdgePair
{
	EdgePair(Triangulation::_Internal::_edge* e1=NULL, Triangulation::_Internal::_edge* e2=NULL, float beta=0, float gamma=0)
	{
		edges[0] = e1;
		edges[1] = e2;
		saliency = 0;
		fitness = 0;
		_temp = 0;
		label = 0;
		count = 0;
		u=NULL; v=NULL; w=NULL;
		if(e1 != NULL && e2 != NULL)
		{
			u = commonVertex(e1, e2);
			v = e1->vertices[0]==u ? e1->vertices[1]: e1->vertices[0];
			w = e2->vertices[0]==u ? e2->vertices[1]: e2->vertices[0];
			fitness = pairFitness(e1, e2, beta, gamma); // * pointFitness(u, v, w, beta, gamma);
			saliency = fitness;
		}
	}
	bool operator <(const EdgePair& pair)
	{
		return saliency < pair.saliency;
	}
	bool operator ==(const EdgePair& pair)
	{
		return (edges[0]==pair.edges[0] && edges[1]==pair.edges[1]) || (edges[0]==pair.edges[1] && edges[1]==pair.edges[0]);
	}
	Triangulation::_Internal::_edge* edges[2];
	Triangulation::_Internal::_vertex* u;
	Triangulation::_Internal::_vertex* v;
	Triangulation::_Internal::_vertex* w;
	float saliency;
	float fitness;
	float _temp;
	int label;
	int count;
};

struct EdgePairSet
{
	EdgePairSet(Triangulation::_Internal::_vertex* u, float beta, float gamma, float thres);
	EdgePairSet(Triangulation::_Internal::_vertex* u0, EdgePair pair)
	{
		u = u0;
		pairs.push_back(pair);
	}
	Triangulation::_Internal::_vertex* u;
	vector<EdgePair> pairs;

	int indexByFitness(Triangulation::_Internal::_edge* e=NULL)
	{
		float max = 0;
		int index = -1;
		for(int i=0; i<pairs.size(); ++i)
		{
			if(e==NULL || pairs[i].edges[0]==e || pairs[i].edges[1]==e)
			{
				if(max < pairs[i].fitness)
				{
					max = pairs[i].fitness;
					index = i;
				}
			}
		}
		return index;
	}
	EdgePair maxBySaliency()
	{
		float maxval = -std::numeric_limits<float>::infinity();
		int index = -1;
		for(int i=0; i<pairs.size(); ++i)
		{
			if(pairs[i].saliency > maxval)
			{
				index = i;
				maxval = pairs[i].saliency;
			}
		}
		if(index>=0)
		{
			return pairs[index];
		}
		else
		{
			return EdgePair(NULL, NULL);
		}
	}
};

EdgePairSet::EdgePairSet(Triangulation::_Internal::_vertex* u0, float beta, float gamma, float thres)
{
	this->u = u0;
	vector<EdgePair> pairs0;
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e1 = u->edges[i];
		for(int j=i+1; j<u->edges.size(); ++j)
		{
			Triangulation::_Internal::_edge* e2 = u->edges[j];
			EdgePair p(e1, e2, beta, gamma);
			pairs0.push_back(p);
		}
	}

	sort(pairs0.begin(), pairs0.end());
	for(int i=0; i<pairs0.size(); ++i)
	{
		int j = pairs0.size()-i-1;
		if(pairs0[j].fitness < thres) break;
		this->pairs.push_back(pairs0[j]);
	}
}

vector<EdgePairSet>
initialize(Triangulation::Triangulator& trmap, float beta, float gamma, float thres)
{
	vector<EdgePairSet> tset;
	for(int i=0; i<trmap.points.size(); ++i)
	{
		tset.push_back(EdgePairSet(trmap.points[i], beta, gamma, thres));
	}
	return tset;
}

#endif /* ___EDGE_PAIR_H___ */