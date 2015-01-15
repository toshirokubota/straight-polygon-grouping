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
#include <TriangulationWithMemento.h>
#include <FragmentInfo.h>
//#include <TriangulationHelper.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
typedef CParticleF graphKey;

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;

struct MovingLineSegment
{
	MovingLineSegment(){}
	MovingLineSegment(CParticleF& p, CParticleF& q)
	{
		this->p = p;
		this->q = q;
	}
	MovingLineSegment(CParticleF& p, CParticleF& q, CParticleF& u)
	{
		this->p = p;
		this->q = q;
		this->perp = u;
	}
	void setPerp(vector<CParticleF>& pnts)
	{
		float x = q.m_X - p.m_X;
		float y = q.m_Y - p.m_Y;
		float len = sqrt(x*x + y*y);
		x /= len;
		y /= len;
		float mx = (p.m_X + q.m_X)/2.0;
		float my = (p.m_Y + q.m_Y)/2.0;
		CParticleF a(mx - y, my + x); //perp to unit
		if(inside(a, pnts))
		{
			perp.m_X = -y;
			perp.m_Y = x;
		}
		else
		{
			perp.m_X = y;
			perp.m_Y = -x;
		}
	}

	//move the line segment in the perpendicular direction
	MovingLineSegment move(float t)
	{
		float x1 = p.m_X + t * perp.m_X;
		float y1 = p.m_Y + t * perp.m_Y;
		float x2 = q.m_X + t * perp.m_X;
		float y2 = q.m_Y + t * perp.m_Y;
		return MovingLineSegment(CParticleF(x1,y1), CParticleF(x2,y2), perp);
	}
	CParticleF p;
	CParticleF q;
	CParticleF perp;
};

bool
coLinear(CParticleF& a, CParticleF& b, CParticleF& c, float precision = 1.0e-3)
{
	float d = Distance2Line(a, c, b);
	return d < precision;
}

vector<CParticleF>
removeDegenerate(vector<CParticleF> points)
{
	while(true)
	{
		vector<CParticleF> result;
		bool bChanged = false;
		for(int i=0; i<points.size(); ++i)
		{
			int i0=(i-1+points.size()) % points.size();
			int i2=(i+1) % points.size();
			if(coLinear(points[i0], points[i], points[i2]))
			{
				bChanged = true;
			}
			else
			{
				result.push_back(points[i]);
			}
		}
		if(bChanged == false)
		{
			return result;
		}
		else
		{
			points = result;
		}
	}
}

CParticleF
intersect(MovingLineSegment& s, MovingLineSegment& t)
{
	pair<float,float> param = _IntersectConvexPolygon::intersect(s.p, s.q, t.p, t.q);
	float a = param.first;
	return CParticleF((1-a)*s.p.m_X + a*s.q.m_X, (1-a)*s.p.m_Y + a*s.q.m_Y);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	//Points
	vector<CParticleF> points;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for(int i=0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(CParticleF(x, y));
		}
	}
	float tmax = 10.0f;
	if(nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(tmax,prhs[1],classMode);
	} 
	float delta = 1.0f;
	if(nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(delta,prhs[2],classMode);
	} 

	points = removeDegenerate(points);

	vector<MovingLineSegment> v;
	for(int i=0; i<points.size(); ++i)
	{
		int i2=(i+1) % points.size();
		v.push_back(MovingLineSegment(points[i], points[i2]));
	}
	for(int i=0; i<points.size(); ++i)
	{
		v[i].setPerp(points);
	}
	int K = (int)(tmax / delta);
	vector<vector<CParticleF>> trace(K+1);
	for(int k=0; k<=K; ++k)
	{
		trace[k] = points;
		for(int i=0; i<v.size(); ++i)
		{
			v[i] = v[i].move(delta);
		}
		for(int i=0; i<v.size(); ++i)
		{
			int i2 = (i+1) % v.size();
			points[i2] = intersect(v[i], v[i2]);
		}
	}
	trace[K] = points;

	if(nlhs >= 1)
	{
		const int dims[] = {points.size(), 2 * (K+1)};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			for(int j=0; j<trace.size(); ++j)
			{
				CParticleF m = trace[j][i];
				SetData2(F, i, 2*j, dims[0], dims[1], m.m_X);
				SetData2(F, i, 2*j+1, dims[0], dims[1], m.m_Y);
			}
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {v.size(), 6};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], v[i].p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], v[i].p.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], v[i].q.m_X);
			SetData2(F, i, 3, dims[0], dims[1], v[i].q.m_Y);
			SetData2(F, i, 4, dims[0], dims[1], v[i].perp.m_X);
			SetData2(F, i, 5, dims[0], dims[1], v[i].perp.m_Y);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	mexUnlock();
}

