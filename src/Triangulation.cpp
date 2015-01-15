#include <Triangulation.h>
#include <algorithm>
#include <cassert>
#include <szParticleF.h>
using namespace std;
#include <TriangulateBourke.h>
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

	/*
	Reconstruct triangulation from indexed set of triangles.
	*/
	void
	Triangulator::_init(const vector<CParticleF>& pnts, const vector<_Internal::_indexed_triangle>& faces)
	{
		for(int i=0; i<pnts.size(); ++i)
		{
			points.push_back(new _Internal::_vertex(pnts[i]));
		}
		for(int i=0; i<faces.size(); ++i)
		{
			_Internal::_vertex* p[3];
			for(int j=0; j<3; ++j)
			{
				p[j] = points[faces[i].index[j]];
			}
			float ccw = CounterClockWise(p[0]->p, p[1]->p, p[2]->p);
			if(ccw < 0)
			{
				//arrange them in clockwise orientation
				swap(p[0], p[2]);
			}
			_Internal::_edge* e[3];
			for(int j=0; j<3; ++j)
			{
				_Internal::_vertex* u = p[j];
				_Internal::_vertex* v = p[(j+1)%3];
				e[j] = new _Internal::_edge(u, v);
				vector<_Internal::_edge*>::iterator it = find_if(edges.begin(), edges.end(), _Internal::_edge_equal(e[j]));
				if(it==edges.end())
				{
					edges.push_back(e[j]);
				}
				else
				{
					delete e[j];
					e[j] = *it;
				}
				u->AddEdge(e[j]);
				v->AddEdge(e[j]);
			}
			_Internal::_triangle* tr = new _Internal::_triangle(e[0], e[1], e[2], p[0], p[1], p[2]);
			this->faces.push_back(tr);
			for(int j=0; j<3; ++j)
			{
				e[j]->AddFace(tr);
			}
		}
		for(int i=0; i<edges.size(); ++i)
		{
			if(edges[i]->faces[0]==NULL && edges[i]->faces[1]==NULL)
			{
				edges[i]->type = _Internal::Isolated;
			}
			else if(edges[i]->faces[0]!=NULL && edges[i]->faces[1]!=NULL)
			{
				edges[i]->type = _Internal::Inside;
			}
			else 
			{
				edges[i]->type = _Internal::Boundary;
				//we do not consider Hole for now...
			}
		}
	}

	Triangulator::Triangulator(const vector<CParticleF>& pnts)
	{
		vector<TriangulateBourke::Triplet> delaunayTriangles = TriangulateBourke::DelaunayTriangulation(pnts);
		vector<_Internal::_indexed_triangle> ifaces;
		for(int i=0; i<delaunayTriangles.size(); i++)
		{
			_Internal::_indexed_triangle itr;
			itr.index[0] = delaunayTriangles[i].index[0];
			itr.index[1] = delaunayTriangles[i].index[1];
			itr.index[2] = delaunayTriangles[i].index[2];
			ifaces.push_back(itr);
		}
		_init(pnts, ifaces);
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
		edges.erase(it);
		deleted_edges.push_back(e);
		//delete e;
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
			}
			_DeleteEdge(e);
			vector<_Internal::_triangle*>::iterator it = find(faces.begin(), faces.end(), tr);
			assert(it != faces.end());
			faces.erase(it);
			deleted_faces.push_back(tr);
			//delete tr;
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
				tr->edges[i]->RemoveFace(tr);
			}
			vector<_Internal::_triangle*>::iterator itf = find(faces.begin(), faces.end(), tr);
			assert(itf != faces.end());
			faces.erase(itf);
			deleted_faces.push_back(tr);
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



