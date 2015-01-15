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
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Kruskal.h>
#include <cassert>
#include <IntersectionConvexPolygons.h>


vector<int>
orderByAngle(vector<CParticleF>& P, vector<int>& J, float x, float y)
{
	vector<pair<float, int>> pairs;
	for (int i = 0; i < J.size(); ++i)
	{
		int j = J[i];
		float a = GetVisualDirection(P[j].m_X, P[j].m_Y, x, y);
		pairs.push_back(pair<float, int>(a, j));
	}
	sort(pairs.begin(), pairs.end());
	vector<int> I;
	for (int i = 0; i < pairs.size(); ++i)
	{
		I.push_back(pairs[i].second);
	}
	return I;
}

float
horizontalDisplacement(CParticleF& a, CParticleF& b, CParticleF& c)
{
	return (c.m_Y - a.m_Y)*(b.m_X - a.m_X) - (c.m_X - a.m_X)*(b.m_Y - a.m_Y);
}

bool
_RosenbergLangridgeCheck(vector<CParticleF>& points, int pi, int qi)
{
	CParticleF p = points[pi];
	CParticleF q = points[qi];
	vector<CParticleF> Plus;
	vector<CParticleF> Minus;
	vector<int> iPlus;
	vector<int> iMinus;
	float dpq = Distance(p, q);
	//if (dpq < eps) return false;
	bool bConnect = true;

	for (int i = 0; i < points.size(); ++i)
	{
		if (i == pi || i == qi) continue;
		CParticleF r = points[i];
		float drp = Distance(r, p);
		float drq = Distance(r, q);
		if (drp <= dpq || drq <= dpq)
		{
			if (horizontalDisplacement(p, q, r) >= 0.0f)
			{
				Plus.push_back(r);
				iPlus.push_back(i);
			}
			else
			{
				Minus.push_back(r);
				iMinus.push_back(i);
			}
		}
	}

	for (int i = 0; i < Plus.size() && bConnect; ++i)
	{
		CParticleF a = Plus[i];
		for (int j = 0; j < Minus.size() && bConnect; ++j)
		{
			CParticleF b = Minus[j];
			float dab = Distance(a, b);
			if (dab < dpq)
			{
				float h1 = horizontalDisplacement(p, b, a);
				if (h1>0)
				{
					float h2 = horizontalDisplacement(b, q, a);
					if (h2>0)
					{
						bConnect = false;
					}
				}
			}
		}
	}
	return bConnect;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: J = RosenbergLangridge(P)");
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

	vector<pair<int, int>> blobs;
	vector<vector<int>> Js;
	for (int i = 0; i < P.size(); ++i)
	{
		vector<int> J;
		for (int j = 0; j < P.size(); ++j)
		{
			if (i == j) continue;
			if (i + 1 == 10 && j + 1 == 12 || i + 1 == 12 && j + 1 == 10)
			{
				i += 0;
			}
			if (_RosenbergLangridgeCheck(P, i, j))
			{
				J.push_back(j);
				blobs.push_back(pair<int, int>(i, j));
			}
		}
		J = orderByAngle(P, J, P[i].m_X, P[i].m_Y);
		Js.push_back(J);
	}
	vector<pair<int, int>> pairs;
	for (int i = 0; i < P.size(); ++i)
	{
		vector<int> J = Js[i];
		for (int j = 0; j < J.size(); ++j)
		{
			int k = J[j];
			int k2 = j == J.size() - 1 ? J[0] : J[j + 1];
			if (find(Js[k].begin(), Js[k].end(), k2) == Js[k].end() &&
				find(Js[k2].begin(), Js[k2].end(), k) == Js[k2].end())
			{
				pairs.push_back(pair<int, int>(i, k));
				//pairs.push_back(pair<int, int>(i, k2));
			}
		}
	}

	if (nlhs >= 1)
	{
		const int dims[] = { pairs.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], pairs[i].first + 1);
			SetData2(F, i, 1, dims[0], dims[1], pairs[i].second + 1);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { blobs.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], blobs[i].first + 1);
			SetData2(F, i, 1, dims[0], dims[1], blobs[i].second + 1);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

