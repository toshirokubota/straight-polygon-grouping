#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>

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
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <Triangulation.h>
#include <FragmentInfo.h>
#include <TriangulationHelper.h>
#include <DistanceMapUtility.h>
#include <ChamferDistance.h>
#include <DotGroupUtility.h>

void
doMedialAxis(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[]);

void
doCentricity(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing Centricity.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: centricity = DotOperations('centricity', dmap)");
		return;
	}
	vector<float> dmap;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		LoadData(dmap, prhs[1], classIdP, ndimP, &dims);
	}
	vector<float> cent(dmap.size(), 0.0f);
	CentricityTransform(cent, dmap, dims);
	if(nlhs >= 1)
	{
		plhs[0] = StoreData(cent, mxSINGLE_CLASS, 2, dims);
	}
}

void
doLocalMaximum(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing Centricity.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: centricity = DotOperations('lmax', dmap)");
		return;
	}
	vector<float> dmap;
	const int* dims;
	{
		mxClassID classIdP;
		int ndimP;
		LoadData(dmap, prhs[1], classIdP, ndimP, &dims);
	}
	float thres = 0.0f;
	if(nrhs>=3)
	{
		mxClassID classMode;
		ReadScalar(thres,prhs[2],classMode);
	}
	vector<CParticleF> lmax = localMaximaPoints(dmap, thres, dims);
	if(nlhs >= 1)
	{
		const int dims[] = {lmax.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], lmax[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], lmax[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
}

void
doTriangulationFromContours(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing triangulation from contours.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: [P F E] = DotOperations('triangulate', contours, interval)");
		return;
	}

	vector<ContourEQW> contours;
	LoadContourEQW(contours, prhs[1]);

	int numSamplePoints = 5;	
	if(nrhs>=3) 
	{
		mxClassID classMode;
		ReadScalar(numSamplePoints,prhs[2],classMode);
	} 

	vector<vector<FragmentInfo>> vinfo = TriangulationPointsFixedInterval(contours, numSamplePoints);
	vector<CParticleF> points;
	rndm(3); //reseed it
	for(int i=0; i<vinfo.size(); ++i)
	{
		for(int j=0; j<vinfo[i].size(); ++j)
		{
			int ci = vinfo[i][j].contourID;
			int pi = vinfo[i][j].pointID;
			CParticleF p(contours[ci].X[pi], contours[ci].Y[pi]);
			if(find(points.begin(), points.end(), p) == points.end())
			{
				points.push_back(p);
			}
		}
	}

	Triangulation::Triangulator trmap(points);

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
}

void
doTriangulationFromPoints(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing triangulation from points.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: [P F E] = DotOperations('delaunay', points)");
		return;
	}

	//Points
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}

	Triangulation::Triangulator trmap(points);

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
}

/*
Compute the Euclidean distance transform of a region enclosed by a sequence of points (SHAPE).
*/
void
doDistanceMap(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing distance transform.\n");
	if(nrhs < 3)
	{
		mexErrMsgTxt("Usage: [D] = DotOperations('dmap', shape, dims)");
		return;
	}
	//Points
	vector<CParticleF> shape;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		shape = vector2particle(P0, dimsP);
	}
	int dims[2];
	{
		vector<int> D;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(D, prhs[2], classIdP, ndimP, &dimsP);
		dims[0] = D[0];
		dims[1] = D[1];
	}

	vector<float> dmap(dims[0]*dims[1], 0.0f);
	updateDistanceMap(dmap, shape, dims);

	if(nlhs >= 1)
	{
		plhs[0] = StoreData(dmap, mxSINGLE_CLASS, 2, dims);
	}
}

/*
Perform connected components of a graph using a disjoint set.
*/
void
doCluster(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing connected component of a graph.\n");
	if(nrhs < 2)
	{
		mexErrMsgTxt("Usage: [P E] = DotOperations('cluster', P, E)");
		return;
	}
	//Points
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}
	//Edges
	vector<int> E;
	const int* dims;
	{
		mxClassID classIdT;
		int ndimT;
		LoadData(E, prhs[2], classIdT, ndimT, &dims);
	}

	vector<Node<int>*> nodes;
	for(int i=0; i<points.size(); ++i) 
	{
		nodes.push_back(makeset(i));
	}
	for(int i=0; i<dims[0]; ++i) 
	{
		int j1 = GetData2(E, i, 0, dims[0], dims[1], 0)-1; //changed to 0 based
		int j2 = GetData2(E, i, 1, dims[0], dims[1], 0)-1; //changed to 0 based
		merge(nodes[j1], nodes[j2]);
	}
	vector<Node<int>*> cl = clusters(nodes);
	vector<int> label(nodes.size());
	for(int i=0; i<label.size(); ++i) 
	{
		int k1 = distance(cl.begin(), find(cl.begin(), cl.end(), findset(nodes[i])));
		SetData2(label, i, 0, label.size(), 1, k1);
	}

	if(nlhs >= 1)
	{
		const int dims[] = {points.size(), 3};
		vector<float> P(dims[0]*dims[1]);
		for(int i=0; i<points.size(); ++i)
		{
			SetData2(P, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(P, i, 1, dims[0], dims[1], points[i].m_Y);
			SetData2(P, i, 2, dims[0], dims[1], (float)label[i]);
		}
		plhs[0] = StoreData(P, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims2[] = {dims[0], 3};
		vector<int> E2(dims2[0]*dims2[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			int j1 = GetData2(E, i, 0, dims[0], dims[1], 0);
			int j2 = GetData2(E, i, 1, dims[0], dims[1], 0);
			SetData2(E2, i, 0, dims2[0], dims2[1], j1);
			SetData2(E2, i, 1, dims2[0], dims2[1], j2);
			SetData2(E2, i, 2, dims2[0], dims2[1], label[j1-1]);
		}
		plhs[1] = StoreData(E2, mxINT32_CLASS, 2, dims2);
	}

	for(int i=0; i<nodes.size(); ++i)
	{
		delete nodes[i];
	}
}

#include <Graph.h>
#include <GraphFactory.h>
#include <Kruskal.h>
GraphFactory<int>* GraphFactory<int>::_instance = 0;
/*
Perform Kruskal's minimum spanning tree.
*/
void
doKruskal(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Doing Kruskal on a graph.\n");
	if (nrhs < 2)
	{
		mexErrMsgTxt("Usage: [P E] = DotOperations('kruskal', P, E, Weight)");
		return;
	}
	//Points
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}
	//Edges
	vector<int> E;
	const int* dims;
	{
		mxClassID classIdT;
		int ndimT;
		LoadData(E, prhs[2], classIdT, ndimT, &dims);
		for (int i = 0; i<E.size(); ++i)
		{
			E[i] -= 1; //work as 0 based index.
		}
	}
	//Weight
	vector<float> W;
	const int* dimsW;
	{
		mxClassID classIdT;
		int ndimT;
		LoadData(W, prhs[3], classIdT, ndimT, &dimsW);
	}

	GraphFactory<int>* factory = GraphFactory<int>::GetInstance();
	vector<Vertex<int>*> vertices;
	vector<Edge<int>*> edges;
	for (int i = 0; i<points.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(i));
	}
	for (int i = 0; i<dims[0]; ++i)
	{
		int k1 = GetData2(E, i, 0, dims[0], dims[1], 0);
		int k2 = GetData2(E, i, 1, dims[0], dims[1], 0);
		float w = GetData2(W, i, 0, dimsW[0], dimsW[1], 0.0f);
		Edge<int>* e1 = factory->makeEdge(vertices[k1], vertices[k2], w);
		Edge<int>* e2 = factory->makeEdge(vertices[k2], vertices[k1], w);
		edges.push_back(e1);
		edges.push_back(e2);
		vertices[k1]->Add(e1);
		vertices[k2]->Add(e2);
	}
	vector<Edge<int>*> mst = Kruskal(edges, vertices);

	if (nlhs >= 1)
	{
		const int dims[] = { points.size(), 3 };
		vector<float> P(dims[0] * dims[1]);
		for (int i = 0; i<points.size(); ++i)
		{
			SetData2(P, i, 0, dims[0], dims[1], points[i].m_X);
			SetData2(P, i, 1, dims[0], dims[1], points[i].m_Y);
		}
		plhs[0] = StoreData(P, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { mst.size(), 2 };
		vector<int> E(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(E, i, 0, dims[0], dims[1], mst[i]->u->key + 1); //change it to one based index
			SetData2(E, i, 1, dims[0], dims[1], mst[i]->v->key + 1);
		}
		plhs[1] = StoreData(E, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { mst.size(), 1 };
		vector<float> E(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(E, i, 0, dims[0], dims[1], mst[i]->w); //change it to one based index
		}
		plhs[2] = StoreData(E, mxSINGLE_CLASS, 2, dims);
	}
	factory->Clean();
}

/*
Calculate the area of a polygon
*/
void
polygonArea(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Calculating the area of a polygon.\n");
	if (nrhs < 2)
	{
		mexErrMsgTxt("Usage: area = DotOperations('area', P)");
		return;
	}
	//Points
	vector<CParticleF> points;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		points = vector2particle(P0, dimsP);
	}

	float area = polygonArea(points);

	if (nlhs >= 1)
	{
		plhs[0] = StoreScalar(area, mxSINGLE_CLASS);
	}
}

vector<CParticleF>
normalizePoints(const vector<CParticleF>& points)
{
	vector<CParticleF> res = points;
	if (points.empty()) return res;

	float maxx = points[0].m_X;
	float minx = maxx;
	float maxy = points[0].m_Y;
	float miny = maxy;
	for (int i = 1; i < points.size(); ++i)
	{
		float x = points[i].m_X;
		float y = points[i].m_Y;
		maxx = Max(maxx, x);
		maxy = Max(maxy, y);
		minx = Min(minx, x);
		miny = Min(miny, y);
	}
	float sx = (maxx - minx) <= 1.0e-5 ? 1.0f : maxx - minx;
	float sy = (maxy - miny) <= 1.0e-5 ? 1.0f : maxy - miny;

	for (int i = 0; i < points.size(); ++i)
	{
		res[i].m_X = (res[i].m_X - minx) / sx;
		res[i].m_Y = (res[i].m_Y - miny) / sy;
	}
	return res;
}

vector<CParticleF>
circularShift(const vector<CParticleF>& points, int shift)
{
	int n = points.size();
	vector<CParticleF> res(n);
	for (int i = 0; i < n; ++i)
	{
		res[i] = points[(i + shift) % n];
	}
	return res;
}

#include <ShapeMatching.h>
/*
Perform shape matching cost calculation.
*/
void
doShapeMatching(int nlhs, mxArray **plhs, int nrhs, const mxArray *prhs[])
{
	printf("Calculating matching cost of two shapes using dynamic programming.\n");
	if (nrhs < 2)
	{
		mexErrMsgTxt("Usage: [cost] = DotOperations('shapematch', P, Q)");
		return;
	}
	//Points 1
	vector<CParticleF> P;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		P = vector2particle(P0, dimsP);
		P = normalizePoints(P);
	}
	//Points 2
	vector<CParticleF> Q;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[2], classIdP, ndimP, &dimsP);
		Q = vector2particle(P0, dimsP);
		Q = normalizePoints(Q);
	}


	if (P.size() > Q.size())
	{
		vector<CParticleF> tmp = P;
		P = Q;
		Q = tmp;
	}
	float mincost = std::numeric_limits<float>::infinity();
	vector<CParticleF> Pmatch = P;
	for (int i = 0; i < P.size(); ++i)
	{
		vector<CParticleF> P2 = circularShift(P, i);
		float c = pathCost(P2, Q);
		if (c < mincost)
		{
			mincost = c;
			Pmatch = P2;
		}
	}

	if (nlhs >= 1)
	{
		plhs[0] = StoreScalar(mincost, mxSINGLE_CLASS);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { Pmatch.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<Pmatch.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], Pmatch[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], Pmatch[i].m_Y);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		const int dims[] = { Q.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<Q.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], Q[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], Q[i].m_Y);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: T = DotOperations(mode)");
		return;
	}

	//Mode
	char strmode[64];
	{
		vector<char> vchars;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(vchars, prhs[0], classIdP, ndimP, &dimsP);
		for(int i=0; i<vchars.size(); ++i) strmode[i] = vchars[i];
		strmode[vchars.size()] = '\0';
	}
	printf("Mode: %s\n", strmode);

	if(stricmp(strmode, "centricity")==0)
	{
		doCentricity(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "lmax")==0)
	{
		doLocalMaximum(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "medaxis")==0)
	{
		doMedialAxis(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "chamfer")==0)
	{
		doChamferDistance(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "triangulate")==0)
	{
		doTriangulationFromContours(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "delaunay")==0)
	{
		doTriangulationFromPoints(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "dmap")==0)
	{
		doDistanceMap(nlhs, plhs, nrhs, prhs);
	}
	else if(stricmp(strmode, "cluster")==0)
	{
		doCluster(nlhs, plhs, nrhs, prhs);
	}
	else if (stricmp(strmode, "kruskal") == 0)
	{
		doKruskal(nlhs, plhs, nrhs, prhs);
	}
	else if (stricmp(strmode, "area") == 0)
	{
		polygonArea(nlhs, plhs, nrhs, prhs);
	}
	else if (stricmp(strmode, "shapematch") == 0)
	{
		doShapeMatching(nlhs, plhs, nrhs, prhs);
	}
	mexUnlock();
}

