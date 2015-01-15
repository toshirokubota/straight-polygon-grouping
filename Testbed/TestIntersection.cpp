#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <Triangulation.h>
#include <FragmentInfo.h>
#include <TriangulationHelper.h>
#include <ConvexityMeasure.h>
#include <intersectionConvexPolygons.h>
using namespace _IntersectConvexPolygon;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	CParticleF a(1.0, 3.0);
	CParticleF b(5.0, 3.0);
	CParticleF c(2.0, 4.0);
	CParticleF d(2.0, 1.0);
	float t = intersect(a, b, c, d);
	float s = intersect(c, d, a, b);
	float u = intersect(a, d, b, c);
	mexUnlock();
}

