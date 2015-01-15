#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>
#include <mex.h>

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szParticleF.h>
#include <FragmentInfo.h>
#include <DistanceMapUtility.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Kruskal.h>

/*
This version allows tracing pixels with the same dvalue.
*/
CParticleF
traceSteepestAscent2(CParticleF p, 
					const vector<float>& dmap, 
					vector<CParticleF>& ridges, 
					const int* dims)
{
	CParticleF q = p;
	vector<CParticleFL> Q(1, p);
	vector<unsigned char> map(dmap.size(), (unsigned char)0);
	while(Q.empty() == false)
	{
		for(int k=0; k<Q.size(); ++k)
		{
			SetData2(map, Q[k].m_X, Q[k].m_Y, dims[0], dims[1], (unsigned char)1);
		}
		vector<CParticleFL> Q2;
		for(int i=0; i<Q.size(); ++i)
		{
			int x = Q[i].m_X;
			int y = Q[i].m_Y;
			float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
			float rmax = 0;
			vector<CParticleFL> nexts;
			for(int k=0; k<NumNeighbors8; ++k)
			{
				int x2 = x + XOffset8[k];
				int y2 = y + YOffset8[k];
				float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
				if(GetData2(map, x2, y2, dims[0], dims[1], (unsigned char)1)==0)
				{
					float len = sqrt((float)XOffset8[k]* XOffset8[k] + YOffset8[k]*YOffset8[k]);
					float rate = (dval2-dval)/len;
					if(rate>= rmax)
					{
						CParticleF q(x2, y2, 0, rate); 
						if(rate > rmax)
						{
							nexts.clear();
						}
						nexts.push_back(q);
						rmax = rate;
					}
				}
			}
			for(int k=0; k<nexts.size(); ++k)
			{
				if(find(Q2.begin(), Q2.end(), nexts[k])==Q2.end())
				{
					Q2.push_back(nexts[k]);
				}
			}
		}
		Q = Q2;
		for(int k=0; k<Q.size(); ++k)
		{
			if(find(ridges.begin(), ridges.end(), Q[k]) != ridges.end())
			{
				q = Q[k];
				Q.clear(); //so that the loop will end...
				break;
			}
		}
	}
	printf("traceSteepestAscent2: running from (%d, %d), reaching (%d, %d)\n", (int)p.m_X, (int)p.m_Y, (int)q.m_X, (int)q.m_Y);
	return q;
}

CParticleF
traceSteepestAscent(CParticleF p, 
					const vector<float>& dmap, 
					vector<CParticleF>& ridges, 
					const int* dims)
{
	int count = 0;
	p.m_Life = 0;
	while(true)
	{
		count++;
		int x = p.m_X;
		int y = p.m_Y;
		float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
		float rmax = 0;
		vector<CParticleF> nexts;
		bool bPeak = true;
		for(int k=0; k<NumNeighbors8; ++k)
		{
			int x2 = x + XOffset8[k];
			int y2 = y + YOffset8[k];
			float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
			float len = sqrt((float)XOffset8[k]* XOffset8[k] + YOffset8[k]*YOffset8[k]);
			float rate = (dval2-dval)/len;
			if(dval2>dval && rate>= rmax)
			{
				CParticleF q(x2, y2, 0, rate); 
				if(rate > rmax)
				{
					nexts.clear();
				}
				nexts.push_back(q);
				rmax = rate;
			}
			if(dval2>=dval) bPeak = false;
		}
		if(nexts.empty())
		{
			if(!bPeak)
			{
				p = traceSteepestAscent2(p, dmap, ridges, dims);
			}
			break; 
		}
		else if(nexts.size()==1) 
		{
			p = nexts[0]; //unambiguous trace found
		}
		else //ambiguous
		{
			p = representativePosition(nexts);
		}
		p.m_Life = count;
		//stop as soon as we reached a ridge point
		if(find(ridges.begin(), ridges.end(), p) != ridges.end())
		{
			break;
		}
	}
	return p;
}

void
doMedialAxis(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing Medial Axis Transform.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: [P mst] = DotOperations('mat', dmap)");
		return;
	}
	vector<float> dmap;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		LoadData(dmap, prhs[1], classIdP, ndimP, &dims);
	}

	vector<float> cent(dmap.size(), 0.0f);
	CentricityTransform(cent, dmap, dims);
	vector<CParticleF> ridges;
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			float c = GetData2(cent, x, y, dims[0], dims[1], 0.0f);
			if(c >= 6.0/8.0)
			{
				ridges.push_back(CParticleF(x, y, 0, c));
			}
		}
	}
	
	const int dimsL[] = {ridges.size(), 3};
	vector<int> links(dims[0] * dims[1]);
	for(int i=0; i<ridges.size(); ++i)
	{
		CParticleF p = traceSteepestAscent(ridges[i], dmap, ridges, dims);
		int j = distance(ridges.begin(), find(ridges.begin(), ridges.end(), p));
		if(j>=ridges.size())
		{
			//Cannot find a leading ridge. Use the closest ridge point as a last resort.
			float dmin = std::numeric_limits<float>::infinity();
			for(int k=0; k<ridges.size(); ++k)
			{
				float d = Distance(p, ridges[k]);
				if(d < dmin)
				{
					dmin = d;
					j = k;
				}
			}
		}
		SetData2(links, i, 0, dimsL[0], dimsL[1], i+1);
		SetData2(links, i, 1, dimsL[0], dimsL[1], j+1);
		SetData2(links, i, 2, dimsL[0], dimsL[1], (int)p.m_Life);
	}

	if(nlhs >= 1)
	{
		//write x, y, centricity value, distance value
		const int dims[] = {ridges.size(), 4};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], ridges[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], ridges[i].m_Y);
			SetData2(F, i, 2, dims[0], dims[1], ridges[i].m_Life);
			float dval = GetData2(dmap, ridges[i].m_X, ridges[i].m_Y, dims[0], dims[1], 0.0f);
			SetData2(F, i, 3, dims[0], dims[1], dval);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		//indices to two ridge points, plus the weight
		plhs[1] = StoreData(links, mxINT32_CLASS, 2, dimsL);
	}
}

