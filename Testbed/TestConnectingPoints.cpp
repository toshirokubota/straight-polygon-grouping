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
#include <TriangulationSaliencyUtility.h>
/*
Sample points along a line connecting the two points.
*/
vector<CParticleF> 
connectPoints(CParticleF& p, CParticleF& q)
{
	vector<CParticleF> line;
	float dx = p.m_X - q.m_X;
	float dy = p.m_Y - q.m_Y;
	float x0, x1, y0, y1;
	if(Abs(dx) > Abs(dy))
	{
		if(p.m_X < q.m_X)
		{
			x0 = p.m_X;
			x1 = q.m_X;
			y0 = p.m_Y;
			y1 = q.m_Y;
		}
		else
		{
			x0 = q.m_X;
			x1 = p.m_X;
			y0 = q.m_Y;
			y1 = p.m_Y;
		}
		float slope = dy / dx;
		for(float x=x0; x<x1; x+=1.0)
		{
			float y = slope * (x - x0) + y0;
			line.push_back(CParticleF(x, y));
		}
	}
	else
	{
		if(p.m_Y < q.m_Y)
		{
			x0 = p.m_X;
			x1 = q.m_X;
			y0 = p.m_Y;
			y1 = q.m_Y;
		}
		else
		{
			x0 = q.m_X;
			x1 = p.m_X;
			y0 = q.m_Y;
			y1 = p.m_Y;
		}
		float slope = dx / dy;
		for(float y=y0; y<y1; y+=1.0)
		{
			float x = slope * (y - y0) + x0;
			line.push_back(CParticleF(x, y));
		}
	}
	line.push_back(CParticleF(x1,y1));
	return line;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [J] = Testbed(points)");
		return;
	}

	//Points
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}
	vector<CParticleF> lines;
	for(int i=0; i<points.size(); ++i)
	{
		int j = (i+1) % points.size();
		vector<CParticleF> line = connectPoints(points[i], points[j]);
		lines.insert(lines.end(), line.begin(), line.end());
	}
	if(nlhs >= 1)
	{
		const int dims[] = {lines.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], lines[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], lines[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

