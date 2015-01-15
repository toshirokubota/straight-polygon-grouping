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
#include <FragmentInfo.h>
#include <DotGroupUtility.h>
#include <TriangulationHelper.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: areas = OverlapArea(P, I, Q, J)");
		return;
	}
	//Points
	vector<CParticleF> P;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		P = vector2particle(P0, dimsP);
	}

	float area = polygonArea(P);
	float len = 0;
	float maxlen = 0;
	vector<float> k(P.size());
	float sumk = 0.0f;
	for (int i = 0; i < P.size(); ++i)
	{
		CParticleF p = P[i];
		CParticleF q = P[(i + 1) % P.size()];
		CParticleF r = P[(i - 1 + P.size()) % P.size()];
		float d = Distance(p, q);
		len += d * d;
		maxlen = Max(d, maxlen);
		float ang = GetVisualAngle(r.m_X, r.m_Y, q.m_X, q.m_Y, p.m_X, p.m_Y);
		float sn = sin(ang/2.0);
		k[i] = 1.0 / (Max(sn*sn, 1.0e-8));
		sumk += k[i];
	}
	//float m = area / len; // sqrt(len);
	vector<float> meas;
	meas.push_back(sqrt(area) / maxlen);

	meas.push_back(sqrt(area) / sqrt(len));
	meas.push_back(meas[1] + 0.001*sumk/k.size());
	if (nlhs >= 1)
	{
		const int dims[] = { meas.size(), 1 };
		plhs[0] = StoreData(meas, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

