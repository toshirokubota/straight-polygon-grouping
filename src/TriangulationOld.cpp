#include <Triangulation.h>
#include <algorithm>
#include <cassert>
#include <szParticleF.h>
using namespace std;
#include <Delaunay.h>
#include <mex.h>

namespace Triangulation
{
	bool isBoundary(const _Internal::_edge* e)
	{
		return e->faces[1] == NULL || e->faces[0] == NULL;
	}
	bool isBoundary(const _Internal::_triangle* t)
	{
		for(int i=0; i<3; ++i)
		{
			if(isBoundary(t->edges[i])) 
				return true;
		}
		return false;
	}

	int tripletSanityCheck(vector<SjaakPriester::Triplet>& triplets, int maxval)
	{
		int violated = 0;
		vector<int> entries;
		for(int i=0; i<triplets.size(); ++i)
		{
			for(int j=0; j<3; ++j)
			{
				int i1 = max(triplets[i].index[j], triplets[i].index[(j+1)%3]);
				int i2 = min(triplets[i].index[j], triplets[i].index[(j+1)%3]);
				int val = i2 * maxval + i1;
				entries.push_back(val);
			}
		}
		sort(entries.begin(), entries.end());
		for(int i=0; i<entries.size(); )
		{
			int count = 0;
			int j = i + 1;
			while(j < entries.size() && entries[j]==entries[i])
			{
				count++;
				j++;
			}
			if(count > 2)
			{
				violated++;
			}
			i = j;
		}
		return violated;
	}

	Triangulator::Triangulator(const vector<CParticleF>& pnts)
	{
		for(int i=0; i<pnts.size(); ++i)
		{
			points.push_back(new _Internal::_vertex(pnts[i]));
		}

		vector<SjaakPriester::Triplet> delaunayTriangles = SjaakPriester::DelaunayTriangulation(pnts);
		/*if(tripletSanityCheck(delaunayTriangles, pnts.size()) > 0)
		{
			int i=0;
		}*/
		for(int i=0; i<delaunayTriangles.size(); i++)
		{
			int xk = delaunayTriangles[i].index[0];
			int yk = delaunayTriangles[i].index[1];
			int zk = delaunayTriangles[i].index[2];
			assert(xk>=0 && xk<pnts.size() && yk>=0 && yk<pnts.size() && zk>=0 && zk<pnts.size());
			CParticleF x = pnts[xk];
			CParticleF y = pnts[yk];
			CParticleF z = pnts[zk];

			float ccw = CounterClockWise(x, y, z);
			if(ccw > 0)
			{
				//arrange them in clockwise orientation
				swap(xk, zk);
				swap(x, z);
			}

			_Internal::_edge* tedges[3];
			tedges[0] = (new _Internal::_edge(points[xk], points[yk]));
			tedges[1] = (new _Internal::_edge(points[yk], points[zk]));
			tedges[2] = (new _Internal::_edge(points[zk], points[xk]));
			_Internal::_triangle* tr = new _Internal::_triangle(0, 0, 0, points[xk], points[yk], points[zk]);
			faces.push_back(tr);
			//add new edges to the edge vector. Make association with a new face and edges.
			for(int j=0; j<3; ++j)
			{
				vector<_Internal::_edge*>::iterator it = find_if(edges.begin(), edges.end(), _Internal::_edge_equal(tedges[j]));
				if(it == edges.end())
				{
					edges.push_back(tedges[j]);
					/*printf("Adding [%2.2f %2.2f]  - [%2.2f %2.2f]\n", 
						max(tedges[j]->vertices[0]->p.m_X, tedges[j]->vertices[0]->p.m_Y),
						min(tedges[j]->vertices[0]->p.m_X, tedges[j]->vertices[0]->p.m_Y),
						max(tedges[j]->vertices[1]->p.m_X, tedges[j]->vertices[1]->p.m_Y),
						min(tedges[j]->vertices[1]->p.m_X, tedges[j]->vertices[1]->p.m_Y));*/
				}
				else
				{
					/*printf("\tNot Adding [%2.2f %2.2f]  - [%2.2f %2.2f]\n", 
						max(tedges[j]->vertices[0]->p.m_X, tedges[j]->vertices[0]->p.m_Y),
						min(tedges[j]->vertices[0]->p.m_X, tedges[j]->vertices[0]->p.m_Y),
						max(tedges[j]->vertices[1]->p.m_X, tedges[j]->vertices[1]->p.m_Y),
						min(tedges[j]->vertices[1]->p.m_X, tedges[j]->vertices[1]->p.m_Y));*/
					delete tedges[j];
					tedges[j] = *it;
				}
				tr->edges[j] = tedges[j];
				tedges[j]->AddFace(tr);
			}
			points[xk]->AddEdge(tedges[0]);
			points[xk]->AddEdge(tedges[2]);
			points[yk]->AddEdge(tedges[0]);
			points[yk]->AddEdge(tedges[1]);
			points[zk]->AddEdge(tedges[1]);
			points[zk]->AddEdge(tedges[2]);
		}
		sort(edges.begin(), edges.end(), _Internal::_edge_less);
		//classify each edge. Initially, they are either Inside or Boundary.
		//After edge-removal, we can have Hole edges when an Inside edge is removed.
		for(int i=0; i<edges.size(); ++i)
		{
			if(isBoundary(edges[i]))
			{
				edges[i]->type = _Internal::Boundary;
			}
			else
			{
				edges[i]->type = _Internal::Inside;
			}
			Node<_Internal::_edge*>* node = makeset(edges[i]);
			vnodes.push_back(node);
			edges[i]->node = node;

			/*printf("%d: [%2.2f %2.2f]  - [%2.2f %2.2f], (%x, %x)\n", i,
				max(edges[i]->vertices[0]->p.m_X, edges[i]->vertices[0]->p.m_Y),
				min(edges[i]->vertices[0]->p.m_X, edges[i]->vertices[0]->p.m_Y),
				max(edges[i]->vertices[1]->p.m_X, edges[i]->vertices[1]->p.m_Y),
				min(edges[i]->vertices[1]->p.m_X, edges[i]->vertices[1]->p.m_Y),
				edges[i]->faces[0], edges[i]->faces[1]);*/
		}
		//make a cluster of boundary edges. Initially, there is only one cluster
		Node<_Internal::_edge*>* root = NULL;
		for(int i=0; i<edges.size(); ++i)
		{
			if(edges[i]->type == _Internal::Boundary)
			{	
				if(root == NULL)
				{
					root = edges[i]->node;
				}
				else
				{
					merge(edges[i]->node, root);
				}
			}
		}
	}

	// Locate a vertex opposite to an edge (E) in a triangle (TR).
	_Internal::_vertex* findOpposite(const _Internal::_triangle* tr, const _Internal::_edge* e)
	{
		if(tr->edges[0]!=e && tr->edges[1]!=e && tr->edges[2]!=e) return NULL;

		for(int i=0; i<3; ++i)
		{
			if(tr->vertices[i]!=e->vertices[0] && tr->vertices[i]!=e->vertices[1])
			{
				return tr->vertices[i];
			}
		}
		return NULL;
	}

	/*
	An edge is deemed safe if removal of it will not make the shape irregular.
	A shape is irregular if a vertex has more than two edges that are on the boundary. 
	*/
	bool isSafe(const _Internal::_edge* e) 
	{
		if(isBoundary(e))
		{
			_Internal::_triangle* tr = e->faces[0] ? e->faces[0] : e->faces[1]; //one has to be non NULL
			if(tr == NULL) return false; //this is a dangling edge. not safe.

			_Internal::_vertex* v = findOpposite(tr, e);
			for(int j=0; j<v->edges.size(); ++j)
			{
				const _Internal::_edge* e2 = v->edges[j];
				if(isBoundary(e2))
				{
					return false;
				}
			}
			return true;
		}
		return false; 
	}
	_Internal::_edge* Triangulator::GetNext() const
	{
		for(int i=edges.size()-1; i>=0; i--)
		{
			if(isSafe(edges[i]))
			{
				return edges[i];
			}
		}
		return NULL;
	}

	int Triangulator::SanityCheck() const
	{
#ifdef _DEBUG
		for(int i=0; i<points.size(); ++i)
		{
			for(int j=0; j<points[i]->edges.size(); ++j)
			{
				if(points[i] != points[i]->edges[j]->vertices[0] && points[i] != points[i]->edges[j]->vertices[1])
				{
					return -1;
				}
			}
		}
		for(int i=0; i<edges.size(); ++i)
		{
			if(edges[i]->faces[0]==NULL && edges[i]->faces[1]==NULL) //a dangling edge should be removed. 
			{
				return (i+1)*1230 + 1;
			}
			if(edges[i]->faces[0])
			{
				if(!(edges[i]->faces[0]->edges[0]==edges[i] || edges[i]->faces[0]->edges[1]==edges[i] || edges[i]->faces[0]->edges[2]==edges[i]))
				{
					return (i+1)*1000 + 2;
				}
			}
			if(edges[i]->faces[1])
			{
				if(!(edges[i]->faces[1]->edges[0]==edges[i] || edges[i]->faces[1]->edges[1]==edges[i] || edges[i]->faces[1]->edges[2]==edges[i]))
				{
					return (i+1)*1000 + 3;
				}
			}
			if(edges[i]->node->key != edges[i])
			{
				return (i+1) * 1000 + 4;
			}
		}
		for(int i=0; i<vnodes.size(); ++i)
		{
			if(vnodes[i]->key->node != vnodes[i])
			{
				return (i+1) * 10000 + 1;
			}
		}
		for(int i=0; i<vnodes.size(); ++i)
		{
			if(find(edges.begin(), edges.end(), vnodes[i]->parent->key)==edges.end())
			{
				return -(i+1) * 10000 + 1;
			}
		}
		for(int i=0; i<faces.size(); ++i)
		{
			for(int j=0; j<3; ++j)
			{
				if(!(faces[i]->edges[j]->faces[0] == faces[i] || faces[i]->edges[j]->faces[1] == faces[i]))
				{
					return -(i+1)*1000 + j;
				}
			}
		}
#endif
		return 0;
	}

	void Triangulator::_DeleteEdge(_Internal::_edge* e)
	{
		vector<_Internal::_edge*>::iterator it = find(edges.begin(), edges.end(), e);
		assert(it != edges.end());
		e->vertices[0]->RemoveEdge(e);
		e->vertices[1]->RemoveEdge(e);
		removeNode(e->node, vnodes);
		edges.erase(it);
		delete e;
	}

	bool Triangulator::RemoveSafeEdge(_Internal::_edge* e)
	{
		if(isSafe(e)) //only a boundary edge can be removed.
		{
			_Internal::_triangle* tr = e->faces[0]; //a triangle to be removed
			for(int i=0; i<3; ++i)
			{
				tr->edges[i]->RemoveFace(tr);
				tr->edges[i]->type == _Internal::Boundary;
				merge(e->node, tr->edges[i]->node);
			}
			_DeleteEdge(e);
			vector<_Internal::_triangle*>::iterator it = find(faces.begin(), faces.end(), tr);
			assert(it != faces.end());
			faces.erase(it);
			delete tr;
		}
		else
		{
			return false;
		}
	}

	bool Triangulator::RemoveIsolatedEdge(_Internal::_edge* e)
	{
		if(e->faces[0]==NULL && e->faces[1]==NULL) 
		{
			_DeleteEdge(e);
			return true;
		}
		else
		{
			return false;
		}
	}

	bool Triangulator::RemoveAnyEdge(_Internal::_edge* e)
	{
		_Internal::_edge* bedge = NULL;
		_Internal::_edge* hedge = NULL;
		if(e->type == _Internal::Inside) e->type = _Internal::Hole; //after removal, a hole will be created.
		for(int k=0; k<2; ++k)
		{
			if(e->faces[k])
			{
				for(int i=0; i<3; ++i)
				{
					if(bedge == NULL && e->faces[k]->edges[i]->type == _Internal::Boundary)
					{
						bedge = e->faces[k]->edges[i];
					}
					if(hedge == NULL && e->faces[k]->edges[i]->type == _Internal::Hole)
					{
						hedge = e->faces[k]->edges[i];
					}
				}
			}
		}
		_Internal::_edge* dominant = e;
		if(bedge) 
		{
			dominant = bedge;
		}
		else if(hedge)
		{
			dominant = hedge;
		}
		for(int k=0; k<2; ++k)
		{
			_Internal::_triangle* tr = e->faces[k]; //a triangle to be removed
			if(tr == NULL) continue;

			//remove this triangle from three edges forming the triangle.
			for(int i=0; i<3; ++i)
			{
				tr->edges[i]->type = dominant->type;
				merge(dominant->node, tr->edges[i]->node);
				tr->edges[i]->RemoveFace(tr);
			}
			vector<_Internal::_triangle*>::iterator itf = find(faces.begin(), faces.end(), tr);
			assert(itf != faces.end());
			faces.erase(itf);
			delete tr;
		}
		for(int i=0; i<edges.size(); ++i)
		{
			if(findset(dominant->node)==findset(edges[i]->node))
			{
				edges[i]->type = dominant->type;
			}
		}
		_DeleteEdge(e);
		return true;
	}

	/*	The result is <0 if the three points are in counter-clock wise.
	It is >0 if they are clock wise.
	It is =0 if they are colinear.
	*/
	int CounterClockWise(const CParticleF& p1, const CParticleF& p2, const CParticleF& p3)
	{
		return (p2.m_X - p1.m_X)*(p3.m_Y - p1.m_Y) - (p2.m_Y - p1.m_Y)*(p3.m_X - p1.m_X);
	}
}



