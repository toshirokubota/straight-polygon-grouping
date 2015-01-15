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

/*
Transform R to the coordinate defined by P and Q where the mid point of P-Q is the origin and P->Q is the x-axis.
*/
CParticleF
transform(CParticleF& r, CParticleF& p, CParticleF& q)
{
	float len = Distance(p, q);
	CParticleF z0((r.m_X - (p.m_X + q.m_X) / 2.0)/len, (r.m_Y - (p.m_Y + q.m_Y) / 2.0)/len); 
	float ang = GetVisualDirection(q.m_X, q.m_Y, p.m_X, p.m_Y);
	float cs = cos(-ang);
	float sn = sin(-ang);
	
	return CParticleF(z0.m_X*cs - z0.m_Y*sn, z0.m_X*sn + z0.m_Y*cs);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: areas = GeometricHashing(P, [scale])");
		return;
	}
	//Points
	map<int,vector<CParticleF>> P;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		if (dimsP[1] < 3)
		{
			mexErrMsgTxt("P has to have at least 4 columns (x, y, model)");
			return;
		}
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], 0.0f);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], 0.0f);
			int model = (int)GetData2(P0, i, 2, dimsP[0], dimsP[1], -1.0f);
			CParticleF p(x, y);
			if (P.find(model) == P.end())
			{
				P[model] = vector<CParticleF>(1, p);
			}
			else
			{
				P[model].push_back(CParticleF(x, y));
			}
		}
	}
	float scale = 1.0;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(scale, prhs[1], classMode);
	}

	float delta = 0.1;
	int count = 0;
	//make a hash table
	map<pair<int,int>, map<int,int>> table;
	for (map<int, vector<CParticleF>>::iterator i = P.begin(); i != P.end(); i++)
	{
		int model = i->first;
		for (vector<CParticleF>::iterator j = i->second.begin(); j != i->second.end(); ++j)
		{
			CParticleF p = *j;
			for (vector<CParticleF>::iterator k = i->second.begin(); k != i->second.end(); ++k)
			{
				CParticleF q = *k;
				if (Distance(p,q) > delta)
				{
					for (vector<CParticleF>::iterator l = i->second.begin(); l != i->second.end(); ++l)
					{
						if (l != j && l != k)
						{
							count++;
							CParticleF r = *l;
							CParticleF z = transform(r, p, q);
							pair<int, int> ref(floor(z.m_X / scale), floor(z.m_Y / scale));
							printf("%d: p=[%2.2f,%2.2f]-q[%2.2f,%2.2f], r[%2.2f,%2.2f] => z[%2.2f,%2.2f], (%d,%d)\n",
								count, p.m_X, p.m_Y, q.m_X, q.m_Y, r.m_X, r.m_Y, z.m_X, z.m_Y, ref.first, ref.second);
							if (table.find(ref) == table.end())
							{
								table[ref] = map<int, int>();
								table[ref][model] = 1;
							}
							else
							{
								if (table[ref].find(model) == table[ref].end())
								{
									table[ref][model] = 1;
								}
								else
								{
									table[ref][model] = table[ref][model] + 1;
								}
							}
						}
					}
				}
			}
		}
	}
	if (nlhs >= 1)
	{
		const int dims[] = { table.size(), 1 };
		vector<vector<int>> vdata(table.size());
		int count = 0;
		for (map<pair<int, int>, map<int, int>>::iterator i = table.begin(); i != table.end(); ++i)
		{
			pair<int,int> ref = i->first;
			vector<int> v;
			v.push_back(ref.first);
			v.push_back(ref.second);
			for (map<int, int>::iterator j = i->second.begin(); j != i->second.end(); ++j)
			{
				v.push_back(j->first);
				v.push_back(j->second);
			}
			vdata[count] = v;
			count++;
		}
		plhs[0] = StoreDataCell(vdata, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

