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
#include <LeastSquaresFitting.h>

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

bool
onSameContourFragment(Triangulation::_Internal::_vertex* u, 
					  Triangulation::_Internal::_vertex* v,
					  vector<FragmentInfo>& vinfo, 
					  vector<CParticleF>& points)
{
	int k1 = distance(points.begin(), find(points.begin(), points.end(), u->p));
	int k2 = distance(points.begin(), find(points.begin(), points.end(), v->p));
	if(vinfo[k1].contourID == vinfo[k2].contourID)
	{
		return true;
	}
	else
	{
		return false;
	}
}

float
geodesicLength(Triangulation::_Internal::_edge* edge, 
				vector<FragmentInfo>& vinfo, 
				vector<CParticleF>& points)
{
	if(onSameContourFragment(edge->vertices[0], edge->vertices[1], vinfo, points)) 
	{
		return 0.0f;
	}
	else
	{
		return edge->Length();
	}
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
			CParticleF p(contours[ci].X[pi], contours[ci].Y[pi]);
			if(find(points.begin(), points.end(), p) == points.end())
			{
				points.push_back(CParticleF(contours[ci].X[pi], contours[ci].Y[pi]));
				vinfo.push_back(vvinfo[i][j]);
			}
		}
	}

	Triangulation::Triangulator trmap(points);
	float slack = 0.01;
	vector<float> vGrad(trmap.edges.size(), 0.0f);
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* edge = trmap.edges[i];
		if(geodesicLength(edge, vinfo, points)<=0)
		{
			continue;
		}
		Triangulation::_Internal::_triangle* tr1 = edge->faces[0];
		Triangulation::_Internal::_triangle* tr2 = edge->faces[1];
		if(tr1==NULL || tr2==NULL)
		{
			continue;
		}
		Triangulation::_Internal::_vertex* op1 = oppositeVertex(tr1, edge);
		Triangulation::_Internal::_vertex* op2 = oppositeVertex(tr2, edge);
		if(onSameContourFragment(op1, edge->vertices[0], vinfo, points) &&
			onSameContourFragment(op2, edge->vertices[0], vinfo, points))
		{
			float a1 = GetVisualAngle(op1->p.m_X, op1->p.m_Y, edge->vertices[0]->p.m_X, edge->vertices[0]->p.m_Y,
				edge->vertices[1]->p.m_X, edge->vertices[1]->p.m_Y);
			float a2 = GetVisualAngle(op2->p.m_X, op2->p.m_Y, edge->vertices[0]->p.m_X, edge->vertices[0]->p.m_Y,
				edge->vertices[1]->p.m_X, edge->vertices[1]->p.m_Y);
			if(a1 + a2 < PI-slack)
			{
				vGrad[i] = 1.0f;
			}
			else if(a1 + a2 > PI+slack)
			{
				vGrad[i] = -1.0f;
			}
		}
		else if (onSameContourFragment(op1, edge->vertices[1], vinfo, points) &&
			onSameContourFragment(op2, edge->vertices[1], vinfo, points))
		{
			float a1 = GetVisualAngle(op1->p.m_X, op1->p.m_Y, edge->vertices[1]->p.m_X, edge->vertices[1]->p.m_Y,
				edge->vertices[0]->p.m_X, edge->vertices[0]->p.m_Y);
			float a2 = GetVisualAngle(op2->p.m_X, op2->p.m_Y, edge->vertices[1]->p.m_X, edge->vertices[1]->p.m_Y,
				edge->vertices[0]->p.m_X, edge->vertices[0]->p.m_Y);
			if(a1 + a2 < PI-slack)
			{
				vGrad[i] = -1.0f;
			}
			else if(a1 + a2 > PI+slack)
			{
				vGrad[i] = 1.0f;
			}
		}
	}
	int numExposed = 0;
	vector<bool> vbExposed(points.size(), false);
	for(int i=0; i<points.size(); ++i)
	{
		for(int j=0; j<trmap.points[i]->edges.size(); ++j)
		{
			if(trmap.points[i]->edges[j]->type == Triangulation::_Internal::Boundary)
			{
				vbExposed[i] = true;
				numExposed++;
				break;
			}
		}
	}
	printf("there are %d exposed points.\n", numExposed);

	const int dimsA[] = {trmap.points.size(), numExposed + trmap.edges.size()};
	vector<double> A(dimsA[0]*dimsA[1], 0.0f);
	vector<double> b(dimsA[1], 0.0f);
	int row = 0;
	double weight = 1.0;
	double weight2 = 10.0;
	for(int i=0; i < points.size(); ++i)
	{
		if(vbExposed[i])
		{
			SetData2(A, i, row, dimsA[0], dimsA[1], weight);
			b[row] = 0;
			row++;
		}
	}
	for(int i=0; i < trmap.edges.size(); ++i)
	{
		int k1 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), trmap.edges[i]->vertices[0]));
		int k2 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), trmap.edges[i]->vertices[1]));
		SetData2(A, k2, row, dimsA[0], dimsA[1], weight2);
		SetData2(A, k1, row, dimsA[0], dimsA[1], -weight2);
		b[row] = vGrad[i] * weight2;
		row++;
	}

	bool bSuccess;
	vector<double> height = LeastSquaresFitting(A, b, points.size(), bSuccess);
	printf("bSuccess = %s, row = %d\n", bSuccess ? "True": "False", row);

	if(nlhs >= 1)
	{
		const int dims[] = {trmap.points.size(), 3};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], points[i].m_Y);
			SetData2(F, i, 2, dims[0], dims[1], (float)height[i]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trmap.faces.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.faces.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[1]->p));
			int j3 = distance(points.begin(), find(points.begin(), points.end(), trmap.faces[i]->vertices[2]->p));
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			SetData2(F, i, 2, dims[0], dims[1], j3+1);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {trmap.edges.size(), 4};
		vector<float> F(dims[0]*dims[1], 0);
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			int j1 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[0]->p));
			int j2 = distance(points.begin(), find(points.begin(), points.end(), trmap.edges[i]->vertices[1]->p));
			SetData2(F, i, 0, dims[0], dims[1], (float)j1+1);
			SetData2(F, i, 1, dims[0], dims[1], (float)j2+1);
			SetData2(F, i, 2, dims[0], dims[1], (float)trmap.edges[i]->type+1);
			SetData2(F, i, 3, dims[0], dims[1], vGrad[i]);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreData(A, mxDOUBLE_CLASS, 2, dimsA);
	}
	if(nlhs >= 5)
	{
		const int dims[] = {dimsA[1], 1};
		plhs[4] = StoreData(b, mxDOUBLE_CLASS, 2, dims);
	}

	mexUnlock();
}
