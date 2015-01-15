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
	//triangle (indices to the points)
	vector<Triangulation::_Internal::_indexed_triangle> T;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
		T = indices2structs(points, T0, dimsT);
	}
	//region
	const int* dimsR;
	vector<int> regions;
	{
		mxClassID classIdT;
		int ndimT;
		LoadData(regions, prhs[2], classIdT, ndimT, &dimsR);
	}
	Triangulation::Triangulator trmap(points, T);
	if(dimsR[0] != trmap.faces.size())
	{
		mexErrMsgTxt("# of regions and # of faces mismatch.");
		return;
	}
	map<Triangulation::_Internal::_edge*,int> strength;
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		strength[trmap.edges[i]] = 0;
	}
	map<Triangulation::_Internal::_triangle*,int> fmap;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		fmap[trmap.faces[i]] = i;
	}
	for(int i=0; i<dimsR[1]; ++i)
	{
		for(int j=0; j<trmap.edges.size(); ++j)
		{
			int c[]={0, 0};
			for(int k=0; k<2; ++k)
			{
				if(trmap.edges[j]->faces[k])
				{
					int id = fmap[trmap.edges[j]->faces[k]];
					if(GetData2(regions, id, i, dimsR[0], dimsR[1], 0) != 0)
					{
						c[k] = 1;
					}
				}
			}
			if(c[0] != c[1])
			{
				strength[trmap.edges[j]]++;
			}
		}
	}
	if(nlhs >= 1)
	{
		const int dims[] = {trmap.edges.size(), 1};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], strength[trmap.edges[i]]);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

