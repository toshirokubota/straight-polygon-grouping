#ifndef ___FLOYD_WARSHALL_H___
#define ___FLOYD_WARSHALL_H___
#include <Graph.h>

AdjacencyPredecessorMatrix
FloydWarshall(const AdjacencyMatrix& A);

vector<int>
FloydWarshall_ReconstructPath(AdjacencyPredecessorMatrix& adj,
							int source, int dest);

#endif /* ___FLOYD_WARSHALL_H___ */