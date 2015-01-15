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

mxArray*
StoreContours(const vector<vector<CParticleF>>& polygons)
{
	const int dims[] = { polygons.size() };
	mxArray* cell = mxCreateCellArray(1, (mwSize*)dims);
	for (int i = 0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size();
		const int dimsC[] = { n, 3 };
		mxArray* ar = mxCreateNumericArray(2, (mwSize*)dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*)mxGetData(ar);
		for (int j = 0; j<n; ++j)
		{
			p[j] = polygons[i][j].m_X;
			p[n + j] = polygons[i][j].m_Y;
			p[2 * n + j] = polygons[i][j].m_Z;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

vector<CParticleF>
reverseShape(vector<CParticleF>& shape)
{
	vector<CParticleF>  rev(shape.size());
	for (int i = 0, j = 0; i < shape.size(); ++i, j = j == 0 ? shape.size() - 1 : j - 1)
	{
		rev[i] = shape[j];
	}
	return rev;
}

CParticleF
triangleCenter(CParticleF& a, CParticleF& b, CParticleF& c)
{
	return CParticleF((a.m_X + b.m_X + c.m_X) / 3.0, (a.m_Y + b.m_Y + c.m_Y) / 3.0, (a.m_Z + b.m_Z + c.m_Z) / 3.0);
}

int
largestOuterTriangle(vector<CParticleF>& vp)
{
	float maxArea = 0;
	int index = -1;
	vector<CParticleF> p(3);
	for (int i = 0; i < vp.size(); ++i)
	{
		p[0] = vp[i == 0 ? vp.size() - 1 : i - 1];
		p[1] = vp[i];
		p[2] = vp[i == vp.size() - 1 ? 0 : i + 1];
		if (ClockWise(p)>0)
		{
			float area = areaTriangle(p[0], p[1], p[2]);
			if (area > maxArea)
			{
				maxArea = area;
				index = i;
			}
		}
	}

	return index;
}

/*
An ear of a polygon is a triangle with two sides being the edges of the polygon and the third one completely inside it.
*/
bool
isEar(int idx, vector<CParticleF>& shape)
{
	CParticleF a = shape[idx == 0 ? shape.size() - 1 : idx - 1];
	CParticleF b = shape[idx];
	CParticleF c = shape[idx == shape.size() - 1 ? 0 : idx + 1];
	CParticleF m((a.m_X + c.m_X) / 2, (a.m_Y + c.m_Y) / 2);
	if (inside(m, shape) == false)
	{
		return false;
	}
	for (int j = 0; j < shape.size(); ++j)
	{
		if (j == idx || j == (idx - 2 + shape.size()) % shape.size() || j == (idx - 1 + shape.size()) % shape.size() || j == (idx + 1) % shape.size()) 
			continue;
		CParticleF p = shape[j];
		CParticleF q = shape[(j + 1) % shape.size()];
		pair<float, float> para = _IntersectConvexPolygon::intersect(a, c, p, q);
		if (para.first > 0 && para.first < 1.0 && para.second > 0 && para.second < 1.0)
		{
			return false;
		}
	}
	return true;
}

/*
This function takes a simple polygon and splitting into triangles using 'ear clipping' method.
*/
vector<vector<CParticleF>>
triangulateInside(vector<CParticleF> shape)
{
	int n = shape.size();
	vector<vector<CParticleF>> triangles;
	while (shape.size() >= 3)
	{
		bool bFound = false;
		for (int i = 0; i < shape.size(); ++i)
		{
			if (shape.size()==3 || isEar(i, shape))
			{
				int idx = i;
				CParticleF a = shape[idx == 0 ? shape.size() - 1 : idx - 1];
				CParticleF b = shape[idx];
				CParticleF c = shape[idx == shape.size() - 1 ? 0 : idx + 1];
				vector<CParticleF> tr(3);
				tr[0] = a; tr[1] = b; tr[2] = c;
				triangles.push_back(tr);
				shape.erase(shape.begin() + idx);
				bFound = true;
				break;
			}
		}
		if (bFound == false)
		{
			/*printf("Incomplete exit: %d triangles out of %d points so far.\n", triangles.size(), n);
			for (int j = 0; j < shape.size(); ++j)
			{
				printf("%f %f\n", shape[j].m_X, shape[j].m_Y);
			}*/
			break;
		}
	}
	return triangles;
}

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
		if (ClockWise(P) < 0)
		{
			P = reverseShape(P);
		}
	}
	//Points
	vector<CParticleF> Q;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		Q = vector2particle(P0, dimsP);
		if (ClockWise(Q) < 0)
		{
			Q = reverseShape(Q);
		}
	}


	vector<vector<CParticleF>> polygon1 = triangulateInside(P);
	vector<vector<CParticleF>> polygon2 = triangulateInside(Q);
	float area1 = 0, area2 = 0, area3 = 0;
	for (int i = 0; i <polygon1.size(); ++i)
	{
		area1 += areaTriangle(polygon1[i][0], polygon1[i][1], polygon1[i][2]);
	}
	for (int i = 0; i <polygon2.size(); ++i)
	{
		area2 += areaTriangle(polygon2[i][0], polygon2[i][1], polygon2[i][2]);
	}
	vector<vector<CParticleF>> overlap;
	for (int i = 0; i < polygon1.size(); ++i)
	{
		for (int j = 0; j < polygon2.size(); ++j)
		{
			vector<CParticleF> cv = _IntersectConvexPolygon::intersectConvexHulls(polygon1[i], polygon2[j]);
			if (cv.size() > 2)
			{
				float a = polygonArea(cv);
				if (a > 0) 
				{
					area3 += a;
					overlap.push_back(cv);
				}
			}
		}
	}

	if(nlhs >= 1)
	{
		const int dims[] = {3, 1};
		vector<float> F(dims[0]*dims[1]);
		F[0] = area1;
		F[1] = area2;
		F[2] = area3;
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreContours(overlap);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreContours(polygon1);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreContours(polygon2);
	}
	mexUnlock();
}

