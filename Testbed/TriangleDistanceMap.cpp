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

vector<float>
calculatePressure(Triangulation::Triangulator& trmap,
				  map<Triangulation::_Internal::_triangle*,float>& dist,
				  map<Triangulation::_Internal::_edge*, float>& border,
				  map<Triangulation::_Internal::_edge*, int>& emap)
{
	vector<pair<float,Triangulation::_Internal::_triangle*>> pairs;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		pairs.push_back(pair<float,Triangulation::_Internal::_triangle*>(dist[trmap.faces[i]], trmap.faces[i]));
	}
	sort(pairs.begin(), pairs.end());
	vector<float> pressure(trmap.edges.size(), 0);
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		if(border[trmap.edges[i]] > 0) pressure[i] = std::numeric_limits<float>::infinity();
	}

	for(int i=0; i<pairs.size(); ++i)
	{
		if(pairs[i].first ==0) continue;

		Triangulation::_Internal::_triangle* tr = pairs[i].second;
		float inflow=0;
		float outflow=0;
		for(int j=0; j<3; ++j)
		{
			Triangulation::_Internal::_edge* ed = tr->edges[j];
			if(border[ed]) continue;
			Triangulation::_Internal::_triangle* tr2 = ed->faces[0]==tr ? ed->faces[1]: ed->faces[0];
			if(tr2)
			{
				if(dist[tr2] < dist[tr])
				{
					inflow += ed->Length();
				}
				else
				{
					outflow += ed->Length();
				}
			}
			else
			{
				outflow += ed->Length();
			}
		}
		for(int j=0; j<3; ++j)
		{
			Triangulation::_Internal::_edge* ed = tr->edges[j];
			if(border[ed]) continue;
			Triangulation::_Internal::_triangle* tr2 = ed->faces[0]==tr ? ed->faces[1]: ed->faces[0];
			if(tr2)
			{
				if(dist[tr2] < dist[tr])
				{
				}
				else
				{
					pressure[emap[ed]] = inflow / outflow;
				}
			}
			else
			{
				pressure[emap[ed]] = inflow / outflow;
			}
		}
	}
	return pressure;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 4 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [J] = Testbed(points, faces, fitness, seed)");
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
	//boundary
	vector<float> boundary;
	{
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(boundary, prhs[2], classIdT, ndimT, &dimsT);
	}
	int seed = -1; //seed triangle
	{
		mxClassID classMode;
		ReadScalar(seed,prhs[3],classMode);
		seed--; //change from MATLAB index.
	} 
	float thres = 10.0;
	if(nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[4],classMode);
	} 
	float height = 10.0;
	if(nrhs >= 6)
	{
		mxClassID classMode;
		ReadScalar(height,prhs[5],classMode);
	} 

	Triangulation::Triangulator trmap(points, T);
	Triangulation::_Internal::_triangle* start = trmap.faces[seed];

	map<Triangulation::_Internal::_vertex*,int> idxmap;
	for(int i=0; i<trmap.points.size(); ++i)
	{
		idxmap[trmap.points[i]] = i;
	}
	map<Triangulation::_Internal::_triangle*,int> fmap;
	map<Triangulation::_Internal::_triangle*,float> fdst;
	map<Triangulation::_Internal::_triangle*,CParticleF> center;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		Triangulation::_Internal::_triangle* tr = trmap.faces[i];
		fmap[tr] = i;
		fdst[tr] = std::numeric_limits<float>::infinity();
		center[trmap.faces[i]] = CParticleF((tr->vertices[0]->p.m_X+tr->vertices[1]->p.m_X+tr->vertices[2]->p.m_X)/3.0,
			(tr->vertices[0]->p.m_Y+tr->vertices[1]->p.m_Y+tr->vertices[2]->p.m_Y)/3.0);
	}
	fdst[start] = 0;
	map<Triangulation::_Internal::_edge*,int> edxmap;
	map<Triangulation::_Internal::_edge*,float> bdrmap;
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		edxmap[trmap.edges[i]] = i;
		bdrmap[trmap.edges[i]] = boundary[i];
	}
	int tic=0;
	while(true)
	{
		tic++;
		bool bChanged = false;
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			Triangulation::_Internal::_triangle* tr = trmap.faces[i];
			CParticleF c1 = center[tr];
			float d = fdst[tr];
			for(int j=0; j<3; ++j)
			{
				Triangulation::_Internal::_edge* ed = tr->edges[j];
				Triangulation::_Internal::_triangle* tr2 = ed->faces[0]==tr ? ed->faces[1]: ed->faces[0];
				if(tr2 != NULL)
				{
					float d2 = fdst[tr2];
					CParticleF c2 = center[tr2];
					float sep = Distance(c1, c2);
					if(bdrmap[ed] <= 0)
					{
						float d3 = d + sep;
						if(d2 > d3)
						{
							fdst[tr2] = d3;
							bChanged = true;
						}
					}
				}
			}
		}
		if(bChanged == false) break;
	}


	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], trmap.points[i]->p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], trmap.points[i]->p.m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trmap.edges.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = idxmap[trmap.edges[i]->vertices[0]];
			int j2 = idxmap[trmap.edges[i]->vertices[1]];
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.faces.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = idxmap[trmap.faces[i]->vertices[0]];
			int j2 = idxmap[trmap.faces[i]->vertices[1]];
			int j3 = idxmap[trmap.faces[i]->vertices[2]];
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		const int dims[] = {trmap.faces.size(), 1};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			if(fdst[trmap.faces[i]] == std::numeric_limits<float>::infinity())
			{
				SetData2(F, i, 0, dims[0], dims[1], -1.0f);
			}
			else
			{
				SetData2(F, i, 0, dims[0], dims[1], fdst[trmap.faces[i]]);
			}
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 5)
	{
		const int dims[] = {trmap.edges.size(), 1};
		vector<float> F(dims[0]*dims[1]);
		vector<float> pressure = calculatePressure(trmap, fdst, bdrmap, edxmap);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], pressure[i]);
		}
		plhs[4] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	mexUnlock();
}

