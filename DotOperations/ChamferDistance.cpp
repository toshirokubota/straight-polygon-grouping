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

void
ChamferDistance(vector<float>& dmap,
				const vector<unsigned char>& fg,
				const int* dims)
{
	for(int i=0; i<fg.size(); ++i)
	{
		if(fg[i]==0) dmap[i] = 0.0f;
		else dmap[i] = std::numeric_limits<float>::infinity();
	}
	int n1 = 3;
	int n2 = 4;
	float nval[] = {n2, n1, n2, n1, n1, n2, n1, n2};
	while(true)
	{
		bool bChanged = false;
		for(int y=0; y<dims[1]; ++y)
		{
			for(int x=0; x<dims[0]; ++x)
			{
				float dval = GetData2(dmap, x, y, dims[0], dims[1], std::numeric_limits<float>::infinity());
				float dmin = dval;
				for(int k=0; k<NumNeighbors8; ++k)
				{
					int x2 = x + XOffset8[k];
					int y2 = y + YOffset8[k];
					float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], std::numeric_limits<float>::infinity());
					if(dval2 + nval[k] < dmin)
					{
						dmin = dval2 + nval[k];
					}
				}
				if(dmin < dval) 
				{
					SetData2(dmap, x, y, dims[0], dims[1], dmin);
					bChanged = true;
				}
			}
		}
		if(bChanged == false) break;
	}
}

void
CenterMaximalBall(vector<unsigned char>& center,
				  const vector<float>& dchamfer,
				  const vector<unsigned char>& fg,
				  const int* dims)
{
	for(int i=0; i<center.size(); ++i)
	{
		center[i] = 0;
	}

	int n1 = 3;
	int n2 = 4;
	float nval[] = {n2, n1, n2, n1, n1, n2, n1, n2};
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(GetData2(fg, x, y, dims[0], dims[1], (unsigned char)1) == 0) continue;

			float dval = GetData2(dchamfer, x, y, dims[0], dims[1], std::numeric_limits<float>::infinity());
			bool bCenter = true;
			for(int k=0; k<NumNeighbors8; ++k)
			{
				int x2 = x + XOffset8[k];
				int y2 = y + YOffset8[k];
				float dval2 = GetData2(dchamfer, x2, y2, dims[0], dims[1], std::numeric_limits<float>::infinity());
				if(dval2 - nval[k] == dval)
				{
					bCenter = false;
					break;
				}
			}
			if(bCenter) 
			{
				SetData2(center, x, y, dims[0], dims[1], (unsigned char)1);
			}
		}
	}
}

void
doChamferDistance(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing Chamfer Distance Transform.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: [D] = DotOperations('chamfer', dmap)");
		return;
	}
	vector<unsigned char> fg;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		LoadData(fg, prhs[1], classIdP, ndimP, &dims);
	}
	vector<float> dmap(fg.size(), 0);
	ChamferDistance(dmap, fg, dims);

	if(nlhs >= 1)
	{
		plhs[0] = StoreData(dmap, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		vector<unsigned char> center(dmap.size(), 0);
		CenterMaximalBall(center, dmap, fg, dims);
		plhs[1] = StoreData(center, mxUINT8_CLASS, 2, dims);
	}
}

