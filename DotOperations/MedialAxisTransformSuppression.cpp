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

pair<float,float>
Gradient(const vector<float>& D,
		 int x, int y,
		 const int* dims)
{
	float vals[8];
	float dval = GetData2(D, x, y, dims[0], dims[1], (float)0);
	for(int k=0; k<NumNeighbors8; ++k)
	{
		int x2 = x + XOffset8[k];
		int y2 = y + YOffset8[k];
		vals[k] = GetData2(D, x2, y2, dims[0], dims[1], (float)dval);
	}
	float gx = ((vals[2]+2*vals[4]+vals[7])-(vals[0]+2*vals[3]+vals[5]))/8.0f;
	float gy = ((vals[5]+2*vals[6]+vals[7])-(vals[0]+2*vals[1]+vals[2]))/8.0f;
	return pair<float,float>(gx, gy);
}

void
FluxTransform(vector<float>& flux, 
			  const vector<float>& dmap, 
			  const int* dims)
{
	vector<float> Gx(dims[0]*dims[1], 0.0f);
	vector<float> Gy(dims[0]*dims[1], 0.0f);

	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			pair<float,float> gr = Gradient(dmap, x, y, dims);
			SetData2(Gx, x, y, dims[0], dims[1], gr.first);
			SetData2(Gy, x, y, dims[0], dims[1], gr.second);
		}
	}

	float nx[]={-1.0f/sqrt(2.0f), 0.0f, 1.0f/sqrt(2.0f), -1.0f, 1.0f, -1.0f/sqrt(2.0f), 0.0f, 1.0f/sqrt(2.0f)};
	float ny[]={-1.0f/sqrt(2.0f), 1.0f, -1.0f/sqrt(2.0f), 0.0f, 0.0f, 1.0f/sqrt(2.0f), 1.0f, 1.0f/sqrt(2.0f)};
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			float sum = 0;
			for(int k=0; k<NumNeighbors8; ++k)
			{
				int x2 = x + XOffset8[k];
				int y2 = y + YOffset8[k];
				float len = sqrt((float)XOffset8[k]*XOffset8[k]+YOffset8[k]*YOffset8[k]);
				float gx = GetData2(Gx, x2, y2, dims[0], dims[1], (float)0);
				float gy = GetData2(Gy, x2, y2, dims[0], dims[1], (float)0);
				sum += (gx*XOffset8[k] + gy*YOffset8[k])/len;
			}
			SetData2(flux, x, y, dims[0], dims[1], sum/NumNeighbors8);
		}
	}
}

bool
isBoundary(const vector<unsigned char>& fg,
		   int x, int y,
		   const int* dims)
{
	if(GetData2(fg, x, y, dims[0], dims[1], (unsigned char)0)==0)
	{
		return false;
	}
	for(int k=0; k<NumNeighbors8; ++k)
	{
		int x2 = x + XOffset8[k];
		int y2 = y + YOffset8[k];
		if(GetData2(fg, x2, y2, dims[0], dims[1], (unsigned char)0))
		{
			return true;
		}
	}
	return false;
}

/*
The pixel (x, y) is simple if it does not break the foreground.
*/
bool
isSimple(const vector<unsigned char>& fg,
		 int x, int y,
		 const int* dims)
{
	if(GetData2(fg, x, y, dims[0], dims[1], (unsigned char)0)==0)
	{
		return false;
	}

	int xoff[]={-1, 0, 1, 1, 1, 0, -1, -1};
	int yoff[]={-1, -1, -1, 0, 1, 1, 1, 0};
	unsigned char vals[8];
	for(int k=0; k<8; ++k)
	{
		int x2 = x + xoff[k];
		int y2 = y + yoff[k];
		vals[k] = GetData2(fg, x2, y2, dims[0], dims[1], (unsigned char)0);
	}

	int count=0; //count 0->1 transitions. there should be only one
	for(int k=0; k<8; ++k)
	{
		if(vals[k]==0 && vals[k==7 ? 0: k+1]==1)
		{
			count++;
		}
	}
	return count==1;
}

bool
isEndPoint(const vector<unsigned char> fg,
		   int x, int y,
		   const int* dims)
{


}

vector<CParticleF>
getSimpleBoundary(const vector<unsigned char> fg,
				  const int* dims)
{
	vector<CParticleF> points;
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(isBoundary(fg, x, y, dims))
			{
				if(isSimple(fg, x, y, dims))
				{
					points.push_back(CParticleF(x, y));
				}
			}
		}
	}
	return points;
}

vector<CParticleF>
medialAxis(vector<unsigned char> fg, //being copied.
		   const vector<float>& flux, 
		   float thres,
		   const int* dims)
{
	vector<CParticleF> Q = getSimpleBoundary(fg, dims);
	for(int i=0; i<Q.size(); ++i)
	{
		Q[i].m_Life = GetData2(flux, Q[i].m_X, Q[i].m_Y, dims[0], dims[1], (float)1.0);
	}
	sort(Q.begin(), Q.end());

	vector<unsigned char> flag(dims[0]*dims[1], (unsigned char)0);
	vector<CParticleF> skeleton;
	while(Q.empty() == false)
	{
		CParticleF q = Q[Q.size()-1];
		int x = q.m_X;
		int y = q.m_Y;
		Q.erase(Q.end()-1);
		SetData2(flag, x, y, dims[0], dims[1], (unsigned char)1);
		if(isSimple(fg, x, y, dims))
		{
			if(q.m_Life < thres)
			{
				skeleton.push_back(q);
			}
			else
			{
				SetData2(fg, x, y, dims[0], dims[1], (unsigned char)0);
				for(int k=0; k<NumNeighbors8; ++k)
				{
					int x2 = x + XOffset8[k];
					int y2 = y + YOffset8[k];
					if(GetData2(flag, x2, y2, dims[0], dims[1], (unsigned char)1)==0)
					{
						if(isSimple(fg, x2, y2, dims))
						{
							float fval = GetData2(flux, x2, y2, dims[0], dims[1], (float)1.0);
							CParticleF q2(x2, y2, 0, fval);
							if(find(Q.begin(), Q.end(), q2)==Q.end())
							{
								Q.push_back(q2);
							}
						}
					}
				}
				sort(Q.begin(), Q.end());
			}
		}
	}
	return skeleton;
}

vector<CParticleF>
collectLocations(const vector<unsigned char>& map, 
				 const int* dims)
{
	vector<CParticleF> locs;
	for(int i=0; i<dims[1]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			if(GetData2(map, j, i, dims[0], dims[1], (unsigned char)0))
			{
				locs.push_back(CParticleF(j, i));
			}
		}
	}
	return locs;
}

#include <NonMaximumSuppression.h>

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
	float thres = -.5f;
	if(nrhs>=3) 
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[2],classMode);
	}

	vector<float> flux(dmap.size(), 0.0f);
	FluxTransform(flux, dmap, dims);
	/*CentricityTransform(flux, dmap, dims);*/

	vector<unsigned char> fg(dmap.size());
	for(int i=0; i<dmap.size(); ++i)
	{
		if(dmap[i]>0) fg[i] = 1;
	}
	//vector<CParticleF> ridges = medialAxis(fg, flux, thres, dims);

	for(int i=0; i<dmap.size(); ++i)
	{
		flux[i] = -flux[i] + 1.0f;
	}
	vector<float> Dx(dims[0]*dims[1], 0.0f);
	vector<float> Dy(dims[0]*dims[1], 0.0f);
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			pair<float,float> gr = Gradient(flux, x, y, dims);
			SetData2(Dx, x, y, dims[0], dims[1], gr.first);
			SetData2(Dy, x, y, dims[0], dims[1], gr.second);
		}
	}
	vector<float> fmat = CannyNonMaximumSuppression(flux, Dx, Dy, dims);
	/*vector<unsigned char> mat(fmat.size(), 0);
	for(int i=0; i<fmat.size(); ++i) mat[i] = fmat[i]>0 ? 1: 0;*/

	vector<unsigned char> mat = thinning(flux, Dx, Dy, dims);
	vector<CParticleF> ridges = collectLocations(mat, dims);

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
		//plhs[0] = StoreData(mat, mxUINT8_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		//indices to two ridge points, plus the weight
		plhs[1] = StoreData(flux, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		//indices to two ridge points, plus the weight
		plhs[2] = StoreData(fg, mxUINT8_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreData(mat, mxUINT8_CLASS, 2, dims);
	}
	if(nlhs >= 5)
	{
		plhs[4] = StoreData(fmat, mxSINGLE_CLASS, 2, dims);
	}
	/*if(nlhs >= 5)
	{
		vector<float> Gx(dims[0]*dims[1], 0.0f);
		vector<float> Gy(dims[0]*dims[1], 0.0f);
		for(int y=0; y<dims[1]; ++y)
		{
			for(int x=0; x<dims[0]; ++x)
			{
				pair<float,float> gr = Gradient(dmap, x, y, dims);
				SetData2(Gx, x, y, dims[0], dims[1], gr.first);
				SetData2(Gy, x, y, dims[0], dims[1], gr.second);
			}
		}
		plhs[3] = StoreData(Gx, mxSINGLE_CLASS, 2, dims);
		plhs[4] = StoreData(Gy, mxSINGLE_CLASS, 2, dims);
	}*/
}

