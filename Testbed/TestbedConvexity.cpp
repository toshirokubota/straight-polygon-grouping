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
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <TriangulationHelper.h>
#include <CircleFitting.h>

Triangulation::_Internal::_vertex*
oppositeVertex(Triangulation::_Internal::_triangle* face,
			   Triangulation::_Internal::_edge* edge)
{
	for(int i=0; i<3; ++i)
	{
		Triangulation::_Internal::_vertex* u = face->vertices[i];
		if(u != edge->vertices[0] && u != edge->vertices[1])
		{
			return u;
		}
	}
	return NULL;
}

/*
Find an edge (u,w) with w!=v.
u and v are from the same contour fragment, and w has to be as well.
The point of w is the closest one to v.
*/
Triangulation::_Internal::_edge* 
adjacentEdge(Triangulation::_Internal::_vertex* u, 
			 Triangulation::_Internal::_vertex* v, 
			 vector<FragmentInfo>& vinfo, 
			 vector<CParticleF>& points)
{
	int ku = distance(points.begin(), find(points.begin(), points.end(), u->p));
	int kv = distance(points.begin(), find(points.begin(), points.end(), v->p));
	int cu = vinfo[ku].contourID;
	int cv = vinfo[kv].contourID;
	if(cu != cv)
	{
		return NULL;
	}
	int pu = vinfo[ku].pointID;
	int pv = vinfo[kv].pointID;
	Triangulation::_Internal::_edge* edge = NULL;
	int pw;
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e2 = u->edges[i];
		Triangulation::_Internal::_vertex* w = e2->vertices[0]==u? e2->vertices[1]: e2->vertices[0];
		int kw = distance(points.begin(), find(points.begin(), points.end(), w->p));
		if(vinfo[kw].contourID == cu)
		{
			if(edge == NULL) 
			{
				edge = e2;
				pw = vinfo[kw].pointID;
			}
			else
			{
				int pw2 = vinfo[kw].pointID;
				if(pu<pv && pv<pw2 && pw2 < pw)
				{
					edge = e2;
					pw = pw2;
				}
				else if(pu>pv && pv>pw2 && pw2>pw)
				{
					edge = e2;
					pw = pw2;
				}
			}
		}
	}
	return edge;
}

/*
Given a triangle (u,v,w), and an edge (w,x) where x!=u and x!=v,
return either (w,x,u) or (w,x,v), whichever is present in the triangulation.
*/
Triangulation::_Internal::_triangle*
adjacentFace2(Triangulation::_Internal::_triangle* face,
			 Triangulation::_Internal::_edge* edge)
{
	for(int i=0; i<2; ++i)
	{
		Triangulation::_Internal::_triangle* tr = edge->faces[i];
		if(tr==NULL) continue;

		int cnt = 0;
		for(int j=0; j<3; ++j)
		{
			if(face->vertices[j]==tr->vertices[0] || face->vertices[j]==tr->vertices[1] || face->vertices[j]==tr->vertices[2])
			{
				cnt++;
			}
		}
		if(cnt == 2)
		{
			return tr;
		}
	}
	return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [] = Testbed(C)");
		return;
	}

	vector<ContourEQW> contours;
	LoadContourEQW(contours, prhs[0]);

	int numSamplePoints = 5;	
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(numSamplePoints,prhs[1],classMode);
	} 

	vector<vector<FragmentInfo>> vvinfo = TriangulationPointsFixedInterval(contours, numSamplePoints);
	vector<CParticleF> points;
	vector<FragmentInfo> vinfo;
	for(int i=0; i<vvinfo.size(); ++i)
	{
		for(int j=0; j<vvinfo[i].size(); ++j)
		{
			int ci = vvinfo[i][j].contourID;
			int pi = vvinfo[i][j].pointID;
			points.push_back(CParticleF(contours[ci].X[pi], contours[ci].Y[pi]));
			vinfo.push_back(vvinfo[i][j]);
		}
	}

	Triangulation::Triangulator trmap(points);
	vector<bool> vConvex(trmap.faces.size(), false);
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		if(i+1==195)
			i+=0;
		Triangulation::_Internal::_triangle* tr = trmap.faces[i];
		int k[3];
		k[0] = distance(points.begin(), find(points.begin(), points.end(), tr->vertices[0]->p));
		k[1] = distance(points.begin(), find(points.begin(), points.end(), tr->vertices[1]->p));
		k[2] = distance(points.begin(), find(points.begin(), points.end(), tr->vertices[2]->p));
		if(vinfo[k[0]].contourID == vinfo[k[1]].contourID && vinfo[k[0]].contourID == vinfo[k[2]].contourID)
		{
			vConvex[i] = true;
		}
		else 
		{
			Triangulation::_Internal::_edge* edge = tr->edges[0];
			for(int j=1; j<3 && vConvex[i]==false; ++j) //try with each edge
			{
				if(edge->Length() > tr->edges[j]->Length())
				{
					edge = tr->edges[j];
				}
			}
			//int k1 = distance(points.begin(), find(points.begin(), points.end(), edge->vertices[0]->p));
			//int k2 = distance(points.begin(), find(points.begin(), points.end(), edge->vertices[1]->p));
			//if(vinfo[k1].contourID == vinfo[k2].contourID && Abs(vinfo[k1].pointID-vinfo[k2].pointID)==1)
			{
				Triangulation::_Internal::_vertex* op1 = oppositeVertex(tr, edge);
				//int k3 = distance(points.begin(), find(points.begin(), points.end(), op1->p));
				for(int k=1; k<2 && vConvex[i]==false; ++k) //try each vertex
				{
					for(int m=0; m<edge->vertices[k]->edges.size() && vConvex[i]==false; ++m) 
					{
						Triangulation::_Internal::_edge* edge2 = edge->vertices[k]->edges[m];
						if(edge2 == edge) continue; //if it is the original edge, skip
						for(int n=0; n<2 && vConvex[i]==false; ++n)
						{
							Triangulation::_Internal::_triangle* tr2 = edge2->faces[n];
							if(tr2 != NULL)
							{
								Triangulation::_Internal::_vertex* op2 = oppositeVertex(tr2, edge2);
								if(op1 == op2)
								{
									vConvex[i] = true;
								}
							}
						}
					}
				}
			}
		}
	}
	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trmap.faces.size(), 4};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[1]->p));
			int j3 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[2]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
			SetData2(F, i, 3, dims[0], dims[1], vConvex[i] ? 1: 0);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.edges.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], (int)trmap.edges[i]->type+1);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	mexUnlock();
}
