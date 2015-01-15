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
#include <szDistanceTransform.h>

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
			if(GetData2(dmap, x, y, dims[0], dims[1], 0.0f)>0)
			{
				for(int k=0; k<NumNeighbors8; ++k)
				{
					int x2 = x + XOffset8[k];
					int y2 = y + YOffset8[k];
					float len = sqrt((float)XOffset8[k]*XOffset8[k]+YOffset8[k]*YOffset8[k]);
					float gx = GetData2(Gx, x2, y2, dims[0], dims[1], (float)0);
					float gy = GetData2(Gy, x2, y2, dims[0], dims[1], (float)0);
					sum += (gx*XOffset8[k] + gy*YOffset8[k])/len;
				}
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
		if(GetData2(fg, x2, y2, dims[0], dims[1], (unsigned char)0)==0)
		{
			return true;
		}
	}
	return false;
}

#include <DisjointSet.h>

/*
return true if the number of 4-neighbor connected components 4-adjacent to (x, y) is 1.
*/
bool
SingleComponent4_4(int values[3][3], int value)
{
	static Node<int>* nodes[3][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	if(nodes[0][0]==NULL)
	{
		for(int i=0; i<3; ++i) 
			for(int j=0; j<3; ++j)
				nodes[i][j] = makeset(3*i+j);
	}
	else
	{
		for(int i=0; i<3; ++i) 
			for(int j=0; j<3; ++j)
				nodes[i][j]->Reset();
	}
	int count = 0;
	for(int i=0; i<3; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			if(values[i][j]==value)
			{
				if(i>0 && values[i-1][j]==value) merge(nodes[i][j], nodes[i-1][j]);
				if(j>0 && values[i][j-1]==value) merge(nodes[i][j], nodes[i][j-1]);
				count++;
			}
		}
	}
	if(count == 0) return false;
	Node<int>* rep = NULL;
	int ioff[]={0, 1, 1, 2};
	int joff[]={1, 0, 2, 1};
	for(int k=0; k<4; ++k)
	{
		int i = ioff[k];
		int j = joff[k];
		if(values[i][j]==value)
		{
			if(rep == NULL) rep = findset(nodes[i][j]);
			else if(rep != findset(nodes[i][j])) return false;
		}
	}
	return true;
}

/*
return true if the number of 8-neighbor connected components 8-adjacent to (x, y) is 1.
*/
bool 
SingleComponent8_8(int values[3][3], int value)
{
	static Node<int>* nodes[3][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	if(nodes[0][0]==NULL)
	{
		for(int i=0; i<3; ++i) 
			for(int j=0; j<3; ++j)
				nodes[i][j] = makeset(3*i+j);
	}
	else
	{
		for(int i=0; i<3; ++i) 
			for(int j=0; j<3; ++j)
				nodes[i][j]->Reset();
	}
	int count = 0;
	for(int i=0; i<3; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			if(values[i][j]==value)
			{
				if(i>0 && j>0 && values[i-1][j-1]==value) merge(nodes[i][j], nodes[i-1][j-1]);
				if(i>0 && values[i-1][j]==value) merge(nodes[i][j], nodes[i-1][j]);
				if(i>0 && j<2 && values[i-1][j+1]==value) merge(nodes[i][j], nodes[i-1][j+1]);
				if(j>0 && values[i][j-1]==value) merge(nodes[i][j], nodes[i][j-1]);
				count++;
			}
		}
	}
	if(count == 0) return false;
	Node<int>* rep = NULL;
	for(int i=0; i<3; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			if(values[i][j]==value)
			{
				if(rep == NULL) rep = findset(nodes[i][j]);
				else if(rep != findset(nodes[i][j])) return false;
			}
		}
	}
	return true;
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
	int vals[3][3];
	for(int i=0; i<3; ++i)
	{
		int y2 = y + i - 1;
		for(int j=0; j<3; ++j)
		{
			int x2 = x + j - 1;
			if(y2>=0 && y2<dims[1] && x2>=0 && x2<dims[0])
			{
				vals[i][j] = GetData2(fg, x2, y2, dims[0], dims[1], (unsigned char)0);
			}
			else
			{
				vals[i][j] = 0;
			}
		}
	}
	vals[1][1]=2;
	bool b1 = SingleComponent8_8(vals, 1);
	//vals[1][1]=0;
	bool b2 = SingleComponent4_4(vals, 0);
	return b1 && b2;
}

bool
isEndPoint(const vector<unsigned char> fg,
		   int x, int y,
		   const int* dims)
{
	if(GetData2(fg, x, y, dims[0], dims[1], (unsigned char)0)==0)
	{
		return false;
	}
	int cnt = 0;
	for(int k=0; k<NumNeighbors8; ++k)
	{
		int x2 = x + XOffset8[k];
		int y2 = y + YOffset8[k];
		if(GetData2(fg, x2, y2, dims[0], dims[1], (unsigned char)0))
		{
			cnt ++;
		}
	}
	return cnt < 2;
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

vector<int>
medialAxis(vector<unsigned char>& fg, 
		   const vector<float>& flux, 
		   float thres,
		   const int* dims)
{
	vector<int> history(dims[0]*dims[1], 0);
	vector<CParticleF> Q = getSimpleBoundary(fg, dims);
	for(int i=0; i<Q.size(); ++i)
	{
		Q[i].m_Life = GetData2(flux, Q[i].m_X, Q[i].m_Y, dims[0], dims[1], (float)1.0);
	}

	vector<CParticleF> ends;
	int iter = 0;
	while(Q.empty() == false)
	{
		iter++;
		sort(Q.begin(), Q.end());
		CParticleF q = Q[Q.size()-1];
		int x = q.m_X;
		int y = q.m_Y;
		Q.erase(Q.end()-1);
		SetData2(history, x, y, dims[0], dims[1], iter);
		if(x+1==38 && y+1==10)
		{
			printf("simple: %d, end: %d, flux: %f, iter=%d\n",
				isSimple(fg, x, y, dims), isEndPoint(fg, x, y, dims), q.m_Life, iter);
			//return history;
		}
		if(isSimple(fg, x, y, dims))
		{
			if(isEndPoint(fg, x, y, dims) && q.m_Life < thres)
			{
			}
			else
			{
				SetData2(fg, x, y, dims[0], dims[1], (unsigned char)0);
				for(int k=0; k<NumNeighbors8; ++k)
				{
					int x2 = x + XOffset8[k];
					int y2 = y + YOffset8[k];
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
			}
		}
	}
	return history;
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
	vector<unsigned char> fg;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		LoadData(fg, prhs[1], classIdP, ndimP, &dims);
	}

	float thres = -.5f;
	if(nrhs>=3) 
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[2],classMode);
	}
	printf("Threshold is %f\n", thres);

	vector<float> dmap(dims[0]*dims[1]);
	vector<unsigned char> bg(dims[0]*dims[1]);
	for(int i=0; i<dims[0]*dims[1]; ++i) 
	{
		bg[i] = fg[i] ? 0: 1;
	}
	DistanceTransformF(dmap, bg, DistanceTransformMode_Euclid, 2, dims);

	vector<float> flux(dmap.size(), 0.0f);
	FluxTransform(flux, dmap, dims);

	vector<int> history = medialAxis(fg, flux, thres, dims);
	vector<CParticleF> ridges;
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(GetData2(fg, x, y, dims[0], dims[1], (unsigned char)0))
			{
				float fval = GetData2(flux, x, y, dims[0], dims[1], 0.0f);
				ridges.push_back(CParticleF(x, y, 0, fval));
			}
		}
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
		plhs[1] = StoreData(fg, mxUINT8_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		plhs[2] = StoreData(flux, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreData(history, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 5)
	{
		plhs[4] = StoreData(dmap, mxSINGLE_CLASS, 2, dims);
	}
}

