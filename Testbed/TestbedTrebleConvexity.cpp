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

Triangulation::_Internal::_edge*
findEdge(Triangulation::Triangulator& trmap,
		 Triangulation::_Internal::_vertex* u, 
		 Triangulation::_Internal::_vertex* v)
{
	for(int i=0; i<u->edges.size(); ++i)
	{
		if(u->edges[i]->vertices[1] == u && u->edges[i]->vertices[0] == v) return u->edges[i];
		else if(u->edges[i]->vertices[0] == u && u->edges[i]->vertices[1] == v) return u->edges[i];
	}
	return NULL;
}

bool convex(Triangulation::_Internal::_vertex* u,
			Triangulation::_Internal::_edge* e1,
			Triangulation::_Internal::_edge* e2,
			Triangulation::_Internal::_edge* e3,
			float slack = 0.05f)
{
	Triangulation::_Internal::_vertex* v;
	Triangulation::_Internal::_vertex* w;
	Triangulation::_Internal::_vertex* x;
	Triangulation::_Internal::_vertex* y;
	if(e1->vertices[0] != e2->vertices[0] && e1->vertices[0] != e2->vertices[1])
	{
		v = e1->vertices[0];
		w = e1->vertices[1];
	}
	else
	{
		v = e1->vertices[1];
		w = e1->vertices[0];
	}
	if(e3->vertices[0] != e2->vertices[0] && e3->vertices[0] != e2->vertices[1])
	{
		y = e3->vertices[0];
		x = e3->vertices[1];
	}
	else
	{
		y = e3->vertices[1];
		x = e3->vertices[0];
	}
	float angle[4];
	angle[0] = GetVisualAngle(v->p.m_X, v->p.m_Y, u->p.m_X, u->p.m_Y, w->p.m_X, w->p.m_Y);
	angle[1] = GetVisualAngle(x->p.m_X, x->p.m_Y, u->p.m_X, u->p.m_Y, w->p.m_X, w->p.m_Y);
	angle[2] = GetVisualAngle(w->p.m_X, w->p.m_Y, u->p.m_X, u->p.m_Y, x->p.m_X, x->p.m_Y);
	angle[3] = GetVisualAngle(y->p.m_X, y->p.m_Y, u->p.m_X, u->p.m_Y, x->p.m_X, x->p.m_Y);
	if(angle[0]+angle[1] < PI-slack && angle[2]+angle[3] < PI-slack)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [] = Testbed(C)");
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

	//trebles
	vector<int> T;
	const int* dims;
	{
		mxClassID classIdT;
		int ndimT;
		LoadData(T, prhs[1], classIdT, ndimT, &dims);
	}

	Triangulation::Triangulator trmap(points);
	vector<vector<int>> trebleEdgeIndices;
	int skipCount = 0;
	for(int i=0; i<dims[0]; ++i)
	{
		int id[4];
		//read point index and convert it to zero based.
		id[0] = (int)GetData2(T, i, 0, dims[0], dims[1], 0) - 1; 
		id[1] = (int)GetData2(T, i, 1, dims[0], dims[1], 0) - 1; 
		id[2] = (int)GetData2(T, i, 2, dims[0], dims[1], 0) - 1; 
		id[3] = (int)GetData2(T, i, 3, dims[0], dims[1], 0) - 1; 
		Triangulation::_Internal::_edge* eid[3];
		if(id[0]<0 || id[1]<0 || id[2]<0 || id[3]<0)
		{
			skipCount++;
			continue;
		}
		eid[0] = findEdge(trmap, trmap.points[id[0]], trmap.points[id[1]]);
		eid[1] = findEdge(trmap, trmap.points[id[1]], trmap.points[id[2]]);
		eid[2] = findEdge(trmap, trmap.points[id[2]], trmap.points[id[3]]);

		if(eid[0]==NULL || eid[1]==NULL || eid[2]==NULL)
		{
			printf("Unmatched edges at row %d", i+1);
			return;
		}
		vector<int> indices;
		indices.push_back(distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), eid[0])));
		indices.push_back(distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), eid[1])));
		indices.push_back(distance(trmap.edges.begin(), find(trmap.edges.begin(), trmap.edges.end(), eid[2])));
		trebleEdgeIndices.push_back(indices);
	}
	printf("Done reading trebles. %d rows were skipped.\n", skipCount);

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

	const int dimsA[] = {trmap.points.size(), numExposed + trebleEdgeIndices.size()};
	vector<double> A(dimsA[0]*dimsA[1], 0.0f);
	vector<double> b(dimsA[1], 0.0f);
	int row = 0;
	double weight = 1.0;
	double weight2 = 1.0;
	for(int i=0; i < points.size(); ++i)
	{
		if(vbExposed[i])
		{
			SetData2(A, i, row, dimsA[0], dimsA[1], weight);
			b[row] = 0;
			row++;
		}
	}

	double unit = 1.0;
	for(int i=0; i < trebleEdgeIndices.size(); ++i)
	{
		Triangulation::_Internal::_edge* edge = trmap.edges[trebleEdgeIndices[i][1]];
		Triangulation::_Internal::_edge* prev = trmap.edges[trebleEdgeIndices[i][0]];
		Triangulation::_Internal::_edge* next = trmap.edges[trebleEdgeIndices[i][2]];
		int k1 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), edge->vertices[0]));
		int k2 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), edge->vertices[1]));
		Triangulation::_Internal::_vertex* op = NULL;
		if(edge->faces[0]) 
		{
			Triangulation::_Internal::_vertex* op0 = oppositeVertex(edge->faces[0], edge);
			if(convex(op0, prev, edge, next))
			{
				op = op0;
			}
		}
		else if(edge->faces[1])
		{
			Triangulation::_Internal::_vertex* op0 = oppositeVertex(edge->faces[1], edge);
			if(convex(op0, prev, edge, next))
			{
				op = op0;
			}
		}
		if(op != NULL)
		{
			int k3 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), op));
			SetData2(A, k1, row, dimsA[0], dimsA[1], weight2/2.0f);
			SetData2(A, k2, row, dimsA[0], dimsA[1], weight2/2.0f);
			SetData2(A, k3, row, dimsA[0], dimsA[1], -weight2);
			b[row] = weight2 * unit;
			row++;
			int k0 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), prev->vertices[0]));
			if(k0==k1)
				k0 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), prev->vertices[1]));
			int k4 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), next->vertices[0]));
			if(k4==k2)
				k4 = distance(trmap.points.begin(), find(trmap.points.begin(), trmap.points.end(), next->vertices[1]));
			/*printf("%d (%d,%d) > %d (%d,%d), %d (%d,%d), %d (%d,%d), %d (%d,%d)\n", 
				k3+1, (int)op->p.m_X, (int)op->p.m_Y, 
				k0+1, (int)points[k0].m_X, (int)points[k0].m_Y, 
				k1+1, (int)points[k1].m_X, (int)points[k1].m_Y, 
				k2+1, (int)points[k2].m_X, (int)points[k2].m_Y, 
				k4+1, (int)points[k4].m_X, (int)points[k4].m_Y);*/
			printf("%d %d %d , %d %d\n", trebleEdgeIndices[i][0]+1, trebleEdgeIndices[i][1]+1, trebleEdgeIndices[i][2]+1,
				(int)op->p.m_X, (int)op->p.m_Y);
		}
		else
		{
			SetData2(A, k1, row, dimsA[0], dimsA[1], weight2);
			SetData2(A, k2, row, dimsA[0], dimsA[1], -weight2);
			b[row] = 0;
			row++;
		}
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
		plhs[1] = StoreData(A, mxDOUBLE_CLASS, 2, dimsA);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {dimsA[1], 1};
		plhs[2] = StoreData(b, mxDOUBLE_CLASS, 2, dims);
	}

	mexUnlock();
}
