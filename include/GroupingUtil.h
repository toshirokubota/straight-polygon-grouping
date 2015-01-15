#ifndef ___GROUPING_UTIL_H___
#define ___GROUPING_UTIL_H___
#include <myDataType.h>
#include <szParticleF.h>
#include <szEdgel.h>
#include <vector>
#include <cmath>
using namespace std;

const int AngleBinSize = 64;

real
Distance(real x, real y, real x2, real y2);

int
findLeftActive(int k, 
			   const vector<EdgelQueue>& elist);

int
findRightActive(int k, 
				const vector<EdgelQueue>& elist);

vector<CParticleF>
OrderEdgePoints(const vector<unsigned char>& I,
				const CParticleF& click,
				const int* dims);

vector<EdgelQueue>
BuildEdgeList(const vector<CParticleF>& vpoints,
			  const CParticleF& click,
			  const int* dims);

vector<EdgelQueue>
BuildEdgeListUnique(const vector<CParticleF>& vpoints,
					const CParticleF& click,
					const int* dims);

int 
numActives(const vector<EdgelQueue>& elist);

vector<Edgel*>
collectActives(const vector<EdgelQueue>& elist);

real
computeFittingError(const vector<Edgel*>& edges,
					const CParticleF& click,
					bool& bFlag);


vector<real> 
computeDistances(const vector<EdgelQueue>&elist, 
				 const CParticleF& click, 
				 const int* dims);

vector<int> 
indexByProximity(const vector<real>& vdistances); 

void
updateOrder(vector<int>& vindex,
			vector<real>& vdistance,
			int k);

#endif /* ___GROUPING_UTIL_H___ */