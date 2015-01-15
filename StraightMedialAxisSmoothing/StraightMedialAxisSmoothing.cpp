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
#include <set>
#include <hash_map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <DisjointSet.h>
//#include <TriangulationWithMemento.h>
//#include <FragmentInfo.h>
//#include <TriangulationHelper.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
typedef CParticleF graphKey;
#include <DeformingPolygonV3.h>
using namespace StraightAxis;
int MovingParticle::_id = 0;
int DeformingPolygon::_id = 0;

struct SnapshotNode
{
	SnapshotNode(int id, int lv=0)
	{
		polygon_id = id;
		level = lv;
	}
	int polygon_id;
	int level;
	vector<CParticleF> outline;
};

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;
GraphFactory<MovingParticle*>* GraphFactory<MovingParticle*>::_instance = NULL;
GraphFactory<SnapshotNode*>* GraphFactory<SnapshotNode*>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;

mxArray*
StoreContours(const vector<vector<CParticleF>>& polygons)
{
	const int dims[] = {polygons.size()};
	mxArray* cell = mxCreateCellArray(1, (mwSize*) dims);
	for(int i=0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size() ;
		const int dimsC[] = {n, 3};
		mxArray* ar = mxCreateNumericArray(2, (mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = polygons[i][j].m_X;
			p[n+j] = polygons[i][j].m_Y;
			p[2*n+j] = polygons[i][j].m_Z;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

vector<pair<int,int>>
indices2pairs(vector<int> T, const int* dims)
{
	int n = dims[0];
	vector<pair<int,int>> pairs(n);
	for(int i=0; i<n; ++i)
	{
		int k1 = GetData2(T, i, 0, dims[0], dims[1], 0);
		int k2 = GetData2(T, i, 1, dims[0], dims[1], 0);
		pairs[i] = pair<int,int>(k1-1, k2-1);
	}
	return pairs;
}

void
printTree(MovingParticle* particle, 
          set<MovingParticle*>& traced,
		  string tab)
{
	if(traced.find(particle) == traced.end())
	{
		traced.insert(particle);
		printf("%s%d(%2.2f, %2.2f)-(%2.2f, %2.2f) %d\n", 
			tab.c_str(), particle->id, particle->p0.m_X, particle->p0.m_Y, particle->p.m_X, particle->p.m_Y, particle->type);
		for(int i=0; i<particle->descendents.size(); ++i)
		{
			printTree(particle->descendents[i], traced, tab + "\t");
		}
	}
}

vector<MovingParticle*>
traceBack(MovingParticle* p)
{
	vector<MovingParticle*> T;
	vector<MovingParticle*> Q(1, p);
	set<MovingParticle*> S;
	S.insert(p);
	while(Q.empty()==false)
	{
		vector<MovingParticle*> Q2;
		for(int i=0; i<Q.size(); ++i)
		{
			if(Q[i]->type == Collide || Q[i]->type==Merge || Q[i]->type==Split)
			{
				for(int j=0; j<Q[i]->descendents.size(); ++j)
				{
					MovingParticle* q = Q[i]->descendents[j];
					if(S.find(q) == S.end())
					{
						Q2.push_back(q);
						S.insert(q);
					}
				}
			}
			else
			{
				T.push_back(Q[i]);
			}
		} 
		Q = Q2;
	}
	return T;
}


vector<CParticleF>
traceBackPolygon(DeformingPolygon& polygon)
{
	vector<CParticleF> trace;
	for(int i=0; i<polygon.particles.size(); ++i)
	{
		MovingParticle* p = polygon.particles[i];
		vector<MovingParticle*> T = traceBack(p);
		vector<pair<float,MovingParticle*>> pairs;
		for(int j=0; j<T.size(); ++j)
		{
			float d = Distance(T[j]->p0, p->prev->p0) - Distance(T[j]->p0, p->next->p0);
			pairs.push_back(pair<float,MovingParticle*>(d, T[j]));
		}
		sort(pairs.begin(), pairs.end());
		for(int j=0; j<pairs.size(); ++j)
		{
			trace.push_back(pairs[j].second->p0);
		}
	}
	return trace;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	//Points
	vector<CParticleF> points;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for(int i=0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(CParticleF(x, y));
			//points.insert(points.begin(), CParticleF(x,y));
		}
	}
	//triangle (indices to the points)
	vector<pair<int,int>> edges;
	{
		vector<int> T0;
		mxClassID classIdT;
		int ndimT;
		const int* dimsT;
		LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
		edges = indices2pairs(T0, dimsT);
	}
	float tthres = 100;
	if(nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(tthres,prhs[2],classMode);
	} 
	int K=1000; //large enough?
	if(nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(K,prhs[3],classMode);
	} 
	float eps = 0.01;
	if(nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(eps,prhs[4],classMode);
		if(eps <=1.0e-6) //too small..
		{
			printf("eps = %f is too small. Reverting it to the default (0.01).\n", eps);
			eps = 0.01;
		}
	} 
	float offset = 0.4f;
	if(nrhs >= 6)
	{
		mxClassID classMode;
		ReadScalar(offset,prhs[5],classMode);
	}
	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	ParticleFactory* pfactory = ParticleFactory::getInstance();

	vector<Vertex<graphKey>*> vertices;
	for(int i=0; i<points.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(points[i]));
	}
	for(int i=0; i<edges.size(); ++i)
	{
		pair<int,int> idx = edges[i];
		//printf("Edge %d - %d\n", idx.first+1, idx.second+1);
		Edge<graphKey>* edge = factory->makeEdge(vertices[idx.first], vertices[idx.second], 1.0f);
		Edge<graphKey>* edge2 = factory->makeEdge(vertices[idx.second], vertices[idx.first], 1.0f);
		vertices[idx.first]->Add(edge);
		vertices[idx.second]->Add(edge2);
	}
	vector<MovingParticle*> trace = traceTree(vertices);
	vector<MovingParticle*> outline = offsetPolygon(trace, offset, 1.0e-5);

	DeformingPolygon dp(0.0f);
	dp.particles = outline;
	vector<DeformingPolygon> polygons;
	polygons.push_back(dp);

	//for creating a tree of snapshots
	GraphFactory<SnapshotNode*>* tfactory = GraphFactory<SnapshotNode*>::GetInstance();
	vector<SnapshotNode*> snapshots;
	snapshots.push_back(new SnapshotNode(dp.id)); //the root
	TreeNode<SnapshotNode*>* root = tfactory->makeTreeNode(snapshots[0]);
	
	vector<vector<CParticleF>> contours;
	float time = 0;
	vector<float> vtimes(1, 0.0f);
	int count = 0;
	while(polygons.empty() == false && time < tthres)
	{
		count++;
		float mint = std::numeric_limits<float>::infinity();
		string eventType = "Unknown";
		MovingParticle* ep = NULL;
		for(int i=0; i<polygons.size(); ++i)
		{
			pair<float,MovingParticle*> e1 = polygons[i].nextEdgeEvent();
			pair<float,MovingParticle*> e2 = polygons[i].nextSplitEvent();
			if(mint > e1.first)
			{
				mint = e1.first;
				eventType = "Edge";
				ep = e1.second;
			}
			if(mint > e2.first)
			{
				mint = e2.first;
				eventType = "Split";
				ep = e2.second;
			}
		}
		time += mint;
		printf("%d: %s %f %f %d %d(%f,%f)\n", count, eventType.c_str(), mint, time, polygons.size(), ep->id, ep->p0.m_X, ep->p0.m_Y);
		vtimes.push_back(time);
		vector<DeformingPolygon> newpolygons;
		for(int i=0; i<polygons.size(); ++i)
		{
			polygons[i].deform(mint);
			vector<DeformingPolygon> polygons2 = fixPolygon(polygons[i], time, eps, false);
			newpolygons.insert(newpolygons.end(), polygons2.begin(), polygons2.end());
			if(polygons2.size() > 1) //branch in the polygon tree
			{
				//snaps.insert(snaps.end(), polygons2.begin(), polygons2.end());
				for(int j=0; j<polygons2.size(); ++j)
				{
					SnapshotNode* snapshot = new SnapshotNode(polygons2[j].id);
					vector<CParticleF> outline = traceBackPolygon(polygons2[j]);
					int ret = ClockWise(outline);
					if(ret < 0)
					{
						snapshot->outline = outline;
					}
					snapshots.push_back(snapshot);
					//look for its parent in the tree.
					TreeNode<SnapshotNode*>* parent = NULL;
					for(int k=0; k<tfactory->vertices.size(); ++k)
					{
						if(tfactory->vertices[k]->key->polygon_id == polygons2[j].parent_id)
						{
							parent = (TreeNode<SnapshotNode*>*) tfactory->vertices[k];
							break;
						}
					}
					if(parent == NULL)
					{
						mexErrMsgTxt("StraightMedialAxisSmoothing: failed to find a parent node for the polygon tree.");
						return;
					}
					TreeNode<SnapshotNode*>* node = tfactory->makeTreeNode(snapshot);
					node->parent = parent;
					parent->children.push_back(node);
				}
			}
		}

		vector<CParticleF> pnts;
		for(int i=0; i<newpolygons.size(); ++i)
		{
			for(int k=0; k<newpolygons[i].particles.size(); ++k)
			{
				CParticleF p = newpolygons[i].particles[k]->p;
				p.m_Z = (float)i; // use Z to separate polygons
				pnts.push_back(p);
			}
		}
		contours.push_back(pnts);

		for(int i=newpolygons.size()-1; i>=0; i--)
		{
			if(newpolygons[i].particles.size() <= 3)
			{
				quickFinish(newpolygons[i].particles, time);
				newpolygons.erase(newpolygons.begin()+i);
			}
		}
		polygons = newpolygons;
	}

	//hierarcyConsistencyCheck();

	if(nlhs >= 1)
	{
		const int dims[] = {outline.size(), 4};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], outline[i]->p0.m_X);
			SetData2(F, i, 1, dims[0], dims[1], outline[i]->p0.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], outline[i]->v[0]);
			SetData2(F, i, 3, dims[0], dims[1], outline[i]->v[1]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {trace.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], trace[i]->p.m_X);
			SetData2(F, i, 1, dims[0], dims[1], trace[i]->p.m_Y);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 3)
	{
		const int dims[] = {pfactory->particles.size(), 8};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			MovingParticle* p = pfactory->particles[i];
			SetData2(F, i, 0, dims[0], dims[1], p->p0.m_X);
			SetData2(F, i, 1, dims[0], dims[1], p->p0.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], p->p.m_X);
			SetData2(F, i, 3, dims[0], dims[1], p->p.m_Y);
			SetData2(F, i, 4, dims[0], dims[1], p->v[0]);
			SetData2(F, i, 5, dims[0], dims[1], p->v[1]);
			SetData2(F, i, 6, dims[0], dims[1], (float)p->id);
			SetData2(F, i, 7, dims[0], dims[1], (float)p->type);
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreContours(contours);
	}
	if(nlhs >= 5)
	{
		vector<vector<CParticleF>> outlines;
		//rearrange in the tree order while collapsing those attached to the empty node to the root
		vector<pair<int,TreeNode<SnapshotNode*>*>> pairs;
		for(int i=0; i<tfactory->vertices.size(); ++i)
		{
			TreeNode<SnapshotNode*>* node = (TreeNode<SnapshotNode*>*) tfactory->vertices[i];
			if(node->key->outline.empty())
			{
				pairs.push_back(pair<int,TreeNode<SnapshotNode*>*>(0, node));
			}
			else
			{
				int level;
				for(int j=0; j<pairs.size(); ++j)
				{
					if(pairs[j].second == node->parent)
					{
						level = pairs[j].first + 1;
					}
				}
				pairs.push_back(pair<int,TreeNode<SnapshotNode*>*>(level, node));
			}
		}
		sort(pairs.begin(), pairs.end());
		for(int i=0; i<pairs.size(); ++i)
		{
			if(pairs[i].second->key->outline.empty()==false && pairs[i].first==1)
			{
				outlines.push_back(pairs[i].second->key->outline);
			}
		}

		plhs[4] = StoreContours(outlines);
	}
	if(nlhs >= 6)
	{
		plhs[5] = StoreStraightAxes(polygons);
	}

	//clean up
	for(int i=0; i<snapshots.size(); ++i)
	{
		delete snapshots[i];
	}
	GraphFactory<graphKey>::GetInstance()->Clean();
	GraphFactory<MovingParticle*>::GetInstance()->Clean();
	GraphFactory<SnapshotNode*>::GetInstance()->Clean();
	pfactory->clean();
	mexUnlock();
}

