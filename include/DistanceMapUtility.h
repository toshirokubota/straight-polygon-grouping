#ifndef ___DISTANCE_MAP_UTILITY_H___
#define ___DISTANCE_MAP_UTILITY_H___
#include <vector>
using namespace std;
#include <szParticleF.h>

vector<CParticleF> localMaximaPoints(const vector<float>& dmap, 
							   float threshold,
							   const int* dims);

vector<CParticleFL> 
clusterLocalMaximaPoints(const vector<CParticleF>& locs,
						 const vector<float>& dmap,
						 float ratio, 
						 const int* dims);

void
CentricityTransform(vector<float>& C, 
					const vector<float>& D,
					const int* dims);

vector<CParticleFL>
traceUpward(CParticleFL p,
			const vector<float>& dmap,
			const int* dims);

vector<CParticleFL>
traceDownward(CParticleFL p,
			const vector<float>& dmap,
			const int* dims);

vector<CParticleF>
traceUpwardRepresentative(CParticleF p,
			const vector<float>& dmap,
			const int* dims);

vector<CParticleF>
traceUpwardPerptoShape(CParticleF p,
			const vector<float>& dmap,
			const CParticleF& adj,
			const int* dims);

vector<CParticleF>
traceUpwardVariableNeighbor(CParticleF p,
			const vector<float>& dmap,
			int maxWidth,
			const int* dims);

vector<CParticleF>
traceUpwardPreferStraight(CParticleF p,
						  const vector<float>& dmap,
						  CParticleF prev,
						  const int* dims);

vector<CParticleFL>
traceDownward(const vector<CParticleFL>& seeds,
			const vector<float>& dmap,
			const int* dims);

#endif /* ___DISTANCE_MAP_UTILITY_H___ */
