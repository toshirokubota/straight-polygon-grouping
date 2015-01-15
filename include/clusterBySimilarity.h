#ifndef ___CLUSTER_BY_SIMILARITY_H___
#define ___CLUSTER_BY_SIMILARITY_H___
#include <vector>
using namespace std;
#include <Triangulation.h>

vector<float>
computeSimilarityMeasures(Triangulation::Triangulator& trmap);

vector<int>
clusterBySimilarity(Triangulation::Triangulator& trmap, float thres);

vector<int>
clusterBySimilarityV2(Triangulation::Triangulator& trmap, 
						float thres, //threshold for adjacent triangle similarity
						float thres2 //threshold for non-adjacent triangle similarity
						);

#endif /* ___CLUSTER_BY_SIMILARITY_H___ */