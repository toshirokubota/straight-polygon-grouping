#include <iostream>
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
#include <set>
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
#include <szDistanceTransform.h>
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>
#include <intersectionConvexPolygons.h>
#include <szConvexHull2D.h>
#include <DeformingPolygonV3.h>
using namespace StraightAxis;

//int MovingParticle::_id = 0;
ParticleFactory* ParticleFactory::_instance = NULL;
int MovingParticle::_id = 0;
int DeformingPolygon::_id = 0;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [J] = Testbed(polygons)");
		return;
	}

	//Points
	vector<DeformingPolygon> polygons;
	LoadStraightAxes(polygons, prhs[0]);
	for(int i=0; i<polygons.size(); ++i)
	{
		int err_code = polygons[i].sanityCheck();
		if(err_code != 0)
		{
			mexErrMsgTxt("Invalid polygon found.");
		}
	}
	if(nlhs >= 1)
	{
		plhs[0] = StoreStraightAxes(polygons);
	}
	
	ParticleFactory* factory = ParticleFactory::getInstance();
	factory->clean();
	mexUnlock();
}

