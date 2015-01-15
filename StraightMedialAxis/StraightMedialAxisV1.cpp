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

CParticleF 
NormalizedDirection(CParticleF& x, CParticleF& o)
{
	float dx = x.m_X - o.m_X;
	float dy = x.m_Y - o.m_Y;
	float len = sqrt(dx*dx + dy*dy);
	if(len > 1.0e-5)
	{
		return CParticleF(dx/len, dy/len);
	}
	else
	{
		return CParticleF(0.0f, 0.0f);
	}
}

struct MovingParticle
{
	MovingParticle(){}
	MovingParticle(vector<CParticleF>& points, int k)
	{
		this->p = points[k];
		int k0 = (k-1+points.size())  % points.size();
		int k2 = (k+1) % points.size();
		CParticleF m((points[k0].m_X+points[k2].m_X)/2, (points[k0].m_Y+points[k2].m_Y)/2);
		this->convex = inside(m, points);
		this->v = _calculateVelocity(p, points[k0], points[k2]);
	}
	CParticleF move(float t)
	{
		p.m_X = p.m_X + t * v.m_X;
		p.m_Y = p.m_Y + t * v.m_Y;
		return p;
	}
	CParticleF _calculateVelocity(CParticleF& o, CParticleF& a, CParticleF& b)
	{
		CParticleF x = NormalizedDirection(a, o);
		CParticleF y = NormalizedDirection(b, o);
		CParticleF z((x.m_X+y.m_X)/2, (x.m_Y+y.m_Y)/2);
		float vx = z.m_X;
		float vy = z.m_Y;
		float len0 = sqrt(vx*vx + vy*vy);
		vx = vx / len0;
		vy = vy / len0;

		float ang = (PI - GetVisualAngle(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y)) / 2.0;
		float len = 1.0f / cos(ang);
		vx = vx * len;
		vy = vy * len;
		if(convex == false)
		{
			vx = -vx;
			vy = -vy;
		}
		return CParticleF(vx, vy);
	}

	CParticleF p;
	CParticleF v;
	bool convex;
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

	vector<MovingParticle> v;
	for(int i=0; i<points.size(); ++i)
	{
		v.push_back(MovingParticle(points, i));
	}

	int K = (int)(tmax / delta);
	vector<vector<CParticleF>> trace(K+1);
	for(int k=0; k<=K; ++k)
	{
		trace[k] = points;
		for(int i=0; i<v.size(); ++i)
		{
			points[i] = v[i].move(delta);
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
		const int dims[] = {v.size(), 4};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], v[i].p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], v[i].p.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], v[i].v.m_X);
			SetData2(F, i, 3, dims[0], dims[1], v[i].v.m_Y);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	mexUnlock();
}

