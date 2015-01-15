#ifndef ___SZ_CONVEXHULL_2D_H___
#define ___SZ_CONVEXHULL_2D_H___
#include <vector>
using namespace std;

#include <szParticleF.h>

vector<CParticleF>
ConvexHull2D(const vector<CParticleF>& vp);

void
PaintInsideHull(vector<unsigned char>& im, 
				const vector<CParticleF>& vch, 
				const int* dims);

/*
Post-processing to remove neraly co-linear points included in a convex hull.
*/
vector<CParticleF>
RemoveColinearPoints(const vector<CParticleF>& hull, float thres=1.0e-5);

#endif /* ___SZ_CONVEXHULL_2D_H___ */