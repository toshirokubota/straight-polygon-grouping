// FloydWarshall.cpp : Defines the entry point for the console application.
//

#include <Graph.h>
#include <limits>
using namespace std;

//set vertices with no edges between to infinity
void
AdjustWeights(AdjacencyMatrix& A)
{
	for(int i=0; i<A.size(); ++i)
	{
		for(int j=0; j<A.size(); ++j)
		{
			if(i == j)
			{
				A.Set(i, j, 0);
			}
			else if(A.Get(i, j) == 0)
			{
				A.Set(i, j, numeric_limits<float>::infinity());
			}
		}
	}
}

//set vertices with no edges between to infinity
void
AdjustPredecessor(AdjacencyPredecessorMatrix& A)
{
	for(int i=0; i<A.size(); ++i)
	{
		for(int j=0; j<A.size(); ++j)
		{
			if(i == j)
			{
				A.SetPi(i, j, -1);
			}
			else if(A.Get(i, j) < numeric_limits<float>::infinity())
			{
				A.SetPi(i, j, i);
			}
		}
	}
}

AdjacencyPredecessorMatrix
FloydWarshall(const AdjacencyMatrix& A)
{
	AdjacencyPredecessorMatrix P = A;
	AdjustPredecessor(P);
	for(int i=0; i<P.size(); ++i)
	{
		for(int j=0; j<P.size(); ++j)
		{
			for(int k=0; k<P.size(); ++k)
			{
				float d0 = P.Get(j,k);
				float d1 = P.Get(j,i);
				float d2 = P.Get(i,k);
				if(d0 > d1 + d2)
				{
					P.Set(j, k, d1 + d2);
					P.SetPi(j, k, P.GetPi(i, k));
					//P.SetPi(i, j, k);
				}
			}
		}
	}

	return P;
}

vector<int>
FloydWarshall_ReconstructPath(AdjacencyPredecessorMatrix& adj,
							int source, int dest)
{
	vector<int> trace;
	if(source < 0 || source >= adj.size() || dest < 0 || dest>= adj.size())
	{
		//return an empty vector.
	}
	else if(adj.GetPi(source, dest)<0)
	{
		trace.push_back(source); //no path from source to dest. Return a vector with source only.
	}
	else //trace the predecessors
	{
		trace.push_back(dest);
		int k = adj.GetPi(source, dest);
		while(k != source)
		{
			trace.insert(trace.begin(), k);
			k = adj.GetPi(source, k);
		}
		trace.insert(trace.begin(), k);
	}
	return trace;
}

