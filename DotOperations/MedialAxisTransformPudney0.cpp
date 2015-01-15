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
#include <ChamferDistance.h>

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
return true if the number of 8-neighbor connected components 4-adjacent to (x, y) is 1.
*/
bool
SingleComponent8_4(int values[3][3], int value)
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
	bool b2 = SingleComponent8_4(vals, 0);
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
	for(int k=0; k<NumNeighbors4; ++k)
	{
		int x2 = x + XOffset4[k];
		int y2 = y + YOffset4[k];
		if(GetData2(fg, x2, y2, dims[0], dims[1], (unsigned char)0))
		{
			cnt ++;
		}
	}
	return cnt < 2;
}

bool
isDeletableA(const vector<unsigned char>& fg,
			const vector<float>& dmap,
			float thres,
		    int x, int y,
		    const int* dims)
{
	return isSimple(fg, x, y, dims) && !isEndPoint(fg, x, y, dims);
}

bool
isDeletableB(const vector<unsigned char>& fg,
			const vector<float>& dmap,
			float thres,
		    int x, int y,
		    const int* dims)
{
	if(isSimple(fg, x, y, dims))
	{
		float dval = GetData2(dmap, x, y, dims[0], dims[1], std::numeric_limits<float>::infinity());
		int n1 = 3;
		int n2 = 4;
		float nval[] = {n2, n1, n2, n1, n1, n2, n1, n2};
		bool bCenter = true;
		for(int k=0; k<NumNeighbors8; ++k)
		{
			int x2 = x + XOffset8[k];
			int y2 = y + YOffset8[k];
			float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], std::numeric_limits<float>::infinity());
			if(dval2 - nval[k] == dval)
			{
				bCenter = false;
				break;
			}
		}
		if(bCenter && dval > thres) 
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}

enum LabelType {Unqueued, Queued, NonDeletable, Unknown};

vector<int>
DistanceOrderedHomotopicThinning(vector<unsigned char>& fg, 
								  const vector<float>& dmap, 
								  bool (*isDeletable)(
										const vector<unsigned char>&, 
										const vector<float>&,
										float,
										int, int, const int*),
								  float thres,
								  const int* dims)
{
	vector<CParticleF> Q;
	vector<LabelType> labels(fg.size());
	vector<int> history(fg.size(), 0);
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(isBoundary(fg, x, y, dims))
			{
				float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
				Q.push_back(CParticleF(x, y, 0, dval));
				SetData2(labels, x, y, dims[0], dims[1], Queued);
			}
			else
			{
				SetData2(labels, x, y, dims[0], dims[1], Unqueued);
			}
		}
	}
	int iter=0;
	while(Q.empty() == false)
	{
		iter++;
		sort(Q.begin(), Q.end());
		CParticleF q = Q[0];
		int x = q.m_X;
		int y = q.m_Y;
		Q.erase(Q.begin());
		SetData2(labels, x, y, dims[0], dims[1], Unqueued);
		if(isDeletable(fg, dmap, thres, x, y, dims))
		{
			SetData2(fg, x, y, dims[0], dims[1], (unsigned char)0);
			SetData2(history, x, y, dims[0], dims[1], iter);
			for(int k=0; k<NumNeighbors8; ++k)
			{
				int x2 = x + XOffset8[k];
				int y2 = y + YOffset8[k];
				{
					if(GetData2(labels, x2, y2, dims[0], dims[1], Unknown) != Queued)
					{
						float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)0);
						Q.push_back(CParticleF(x2, y2, 0, dval2));
						SetData2(labels, x2, y2, dims[0], dims[1], Queued);
					}
				}
			}
		}
		else
		{
			SetData2(labels, x, y, dims[0], dims[1], NonDeletable);
		}
	}
	return history;
}

void
doMedialAxis(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing Medial Axis Transform Pudney.\n");
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
	float thres = 5.0f;
	if(nrhs>=3) 
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[2],classMode);
	}
	vector<float> dmap(fg.size());
	ChamferDistance(dmap, fg, dims);

	vector<int> labels = DistanceOrderedHomotopicThinning(fg, dmap, *isDeletableB, thres, dims);
	labels = DistanceOrderedHomotopicThinning(fg, dmap, *isDeletableA, thres, dims);

	vector<CParticleF> ridges;
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(GetData2(fg, x, y, dims[0], dims[1], (unsigned char) 0)>0)
			{
				ridges.push_back(CParticleF(x, y, 0, GetData2(dmap, x, y, dims[0], dims[1], (float)0)));
			}
		}
	}

	if(nlhs >= 1)
	{
		//write x, y, centricity value, distance value
		const int dims[] = {ridges.size(), 3};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], ridges[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], ridges[i].m_Y);
			SetData2(F, i, 2, dims[0], dims[1], ridges[i].m_Life);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StoreData(fg, mxUINT8_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		plhs[2] = StoreData(dmap, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreData(labels, mxINT32_CLASS, 2, dims);
	}
}

