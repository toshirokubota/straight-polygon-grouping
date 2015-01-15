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
#include <map>
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
#include <ConvexityMeasure.h>
#include <intersectionConvexPolygons.h>
using namespace _IntersectConvexPolygon;

vector<int>
separateTree(vector<CParticleF>& points,
			 vector<int>& tree,
			 int removeId,
			 const int* dims)
{
	vector<Node<int>*> nodes;
	for(int i=0; i<points.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for(int i=0; i<dims[0]; ++i)
	{
		if(i != removeId)
		{
			int j1 = GetData2(tree, i, 0, dims[0], dims[1], -1);
			int j2 = GetData2(tree, i, 1, dims[0], dims[1], -1);
			merge(nodes[j1], nodes[j2]);
		}
	}
	vector<Node<int>*> cl = clusters(nodes);
	vector<int> labels(points.size(), 0);
	for(int i=0; i<nodes.size(); ++i)
	{
		int k = distance(cl.begin(), find(cl.begin(), cl.end(), findset(nodes[i])));
		labels[i] = k + 1;
	}
	for(int i=0; i<nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}

vector<CParticleF>
convexHullFromTriangles(vector<Triangulation::_Internal::_triangle*>& faces)
{
	vector<CParticleF> points;
	for(int i=0; i<faces.size(); ++i)
	{
		for(int j=0; j<3; ++j)
		{
			if(find(points.begin(), points.end(), faces[i]->vertices[j]->p)==points.end())
			{
				points.push_back(faces[i]->vertices[j]->p);
			}
		}
	}
	vector<CParticleF> hull = ConvexHull2D(points);
	return hull;
}

float
AreaTriangulated(vector<Triangulation::_Internal::_triangle*>& left)
{
	float sum = 0;
	for(int i=0; i<left.size(); ++i)
	{
		sum += areaTriangle(left[i]->vertices[0]->p, left[i]->vertices[1]->p, left[i]->vertices[2]->p);
	}
	return sum;
}

float
AreaConvexHull(vector<CParticleF>& hull)
{
	if(hull.size() < 3) return 0.0f;

	float sum = 0;
	for(int i=1; i<hull.size()-1; ++i)
	{
		sum += areaTriangle(hull[0], hull[i], hull[i+1]);
	}
	return sum;
}

float 
separation(vector<Triangulation::_Internal::_triangle*>& left,
		   vector<Triangulation::_Internal::_triangle*>& right,
		   vector<Triangulation::_Internal::_triangle*>& outside)
{
	vector<CParticleF> lhull = convexHullFromTriangles(left);
	vector<CParticleF> rhull = convexHullFromTriangles(right);
	//vector<CParticleF> ohull = convexHullFromTriangles(outside);
	float areaTL = AreaTriangulated(left);
	float areaCL = AreaConvexHull(lhull);
	float areaTR = AreaTriangulated(right);
	float areaCR = AreaConvexHull(rhull);
	//float areaTO = AreaTriangulated(outside);
	//float areaCO = AreaConvexHull(ohull);

	//vector<CParticleF> cap = intersectConvexHulls(lhull, rhull);
	float minArea = (areaTL+areaTR) / 5.0; //penalize a small area
	float measure = areaTL/Max(areaCL, minArea) * areaTR/Max(areaCR, minArea);// * areaTO/Max(areaCO, minArea);
	return measure;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 3 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: TestBed(P, tree, id)");
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
	const int* dimsF;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsF);
		T = indices2structs(points, T0, dimsF);
	}
	//spanning tree
	vector<int> tree;
	const int* dimsT;
	{
		mxClassID classIdT;
		int ndimT;
		LoadData(tree, prhs[2], classIdT, ndimT, &dimsT);
		for(int i=0; i<tree.size(); ++i) tree[i]--; //one index to zero index
	}

	int removeEdge = -1.0;	
	{
		mxClassID classMode;
		ReadScalar(removeEdge,prhs[3],classMode);
		removeEdge-=1; //one index to zero index
	} 
	Triangulation::Triangulator trmap(points, T);
	map<Triangulation::_Internal::_vertex*,int> lut; //to facilitate point-index lookup
	for(int i=0; i<trmap.points.size(); ++i)
	{
		lut[trmap.points[i]] = i;
	}

	vector<int> label = separateTree(points, tree, removeEdge, dimsT);
	vector<Triangulation::_Internal::_triangle*> outside;
	vector<Triangulation::_Internal::_triangle*> left;
	vector<Triangulation::_Internal::_triangle*> right;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		int k1 = lut[trmap.faces[i]->vertices[0]];
		int k2 = lut[trmap.faces[i]->vertices[1]];
		int k3 = lut[trmap.faces[i]->vertices[2]];
		if(label[k1] == 1 && label[k2] == 1 && label[k3]== 1)
		{
			left.push_back(trmap.faces[i]);
		}
		else if(label[k1] == 2 && label[k2] == 2 && label[k3]== 2)
		{
			right.push_back(trmap.faces[i]);
		}
		else
		{
			outside.push_back(trmap.faces[i]);
		}
	}

	float measure = separation(left, right, outside);
	printf("%d %d %f\n", removeEdge+1, outside.size(), measure);

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
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 3)
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
			if(find(left.begin(), left.end(), trmap.faces[i])<left.end())
			{
				SetData2(F, i, 3, dims[0], dims[1], 1);
			}
			else if(find(right.begin(), right.end(), trmap.faces[i])<right.end())
			{
				SetData2(F, i, 3, dims[0], dims[1], 2);
			}
			else
			{
				SetData2(F, i, 3, dims[0], dims[1], 0);
			}
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		const int dims[] = {dimsT[0], 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<dims[0]; ++i)
		{
			int j1 = GetData2(tree, i, 0, dimsT[0], dimsT[1], -1);
			int j2 = GetData2(tree, i, 1, dimsT[0], dimsT[1], -1);
			SetData2(F, i, 0, dims[0], dims[1], j1+1);
			SetData2(F, i, 1, dims[0], dims[1], j2+1);
			if(i == removeEdge)
			{
				SetData2(F, i, 2, dims[0], dims[1], (int)0);
			}
			else 
			{
				SetData2(F, i, 2, dims[0], dims[1], label[j1]);
			}
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if(nlhs >= 5)
	{
		plhs[4] = StoreScalar(measure, mxSINGLE_CLASS);
	}
	/*if(nlhs >= 5)
	{
		const int dims[] = {lhull.size()+rhull.size()+cap.size(), 3};
		vector<float> F(dims[0]*dims[1], 0);
		int j = 0;
		for(int i=0; i<lhull.size(); ++i, ++j)
		{
			SetData2(F, j, 0, dims[0], dims[1], lhull[i].m_X);
			SetData2(F, j, 1, dims[0], dims[1], lhull[i].m_Y);
			SetData2(F, j, 2, dims[0], dims[1], 1.0f);
		}
		for(int i=0; i<rhull.size(); ++i, ++j)
		{
			SetData2(F, j, 0, dims[0], dims[1], rhull[i].m_X);
			SetData2(F, j, 1, dims[0], dims[1], rhull[i].m_Y);
			SetData2(F, j, 2, dims[0], dims[1], 2.0f);
		}
		for(int i=0; i<cap.size(); ++i, ++j)
		{
			SetData2(F, j, 0, dims[0], dims[1], cap[i].m_X);
			SetData2(F, j, 1, dims[0], dims[1], cap[i].m_Y);
			SetData2(F, j, 2, dims[0], dims[1], 3.0f);
		}
		plhs[4] = StoreData(F, mxINT32_CLASS, 2, dims);
	}*/

	mexUnlock();
}

