#include <TriangulationWithMemento.h>
#include <algorithm>
#include <cassert>
//#include <szParticleF.h>
using namespace std;
#include <TriangulateBourke.h>
#include <mex.h>

namespace Triangulation
{
	namespace _Internal
	{
		bool isBoundary(const _edge* e)
		{
			return e->faces[0]==NULL || e->faces[1]==NULL || isActive(e->faces[1])==false || isActive(e->faces[0]) == false;
		}
		bool isBoundary(const _triangle* t)
		{
			for(int i=0; i<3; ++i)
			{
				if(isBoundary(t->edges[i])) 
					return true;
			}
			return false;
		}

		bool isActive(const _edge* e) {return e->bActive;}
		bool isActive(const _triangle* t) {return t!=NULL && t->bActive;}
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
			if(ccw > 0)
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


	void Triangulator::_DeleteEdge(_Internal::_edge* e)
	{
		e->bActive = false;
		e->vertices[0]->RemoveEdge(e);
		e->vertices[1]->RemoveEdge(e);
	}

	bool Triangulator::RemoveEdge(_Internal::_edge* e)
	{
		_Internal::_Memento m;
		m.edge = e;
		m.face = NULL;
		if(isBoundary(e) == false) return false;
		if(_Internal::isActive(e->faces[0]) == false && _Internal::isActive(e->faces[1])) 
		{
			e->faces[1]->bActive = false;
			m.face = e->faces[1];
		}
		else if(_Internal::isActive(e->faces[1]) == false && _Internal::isActive(e->faces[0]))
		{
			e->faces[0]->bActive = false;
			m.face = e->faces[0];
		}
		_DeleteEdge(e);
		states.push_back(m);
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



