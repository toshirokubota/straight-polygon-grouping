#ifndef ___DOT_GROUP_UTILITY_H___
#define ___DOT_GROUP_UTILITY_H___

#include <vector>
#include <szParticleF.h>
#include <Triangulation.h>

bool
inside(const CParticleF& p, const vector<CParticleF>& pnts);

void
updateDistanceMap(vector<float>& dmap, 
				  const vector<CParticleF>& shape, 
				  const int* dims);

#endif /* ___DOT_GROUP_UTILITY_H___ */