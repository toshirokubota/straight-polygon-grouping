#ifndef ___TRIANGULATION_WITH_MEMENTO_H___
#define ___TRIANGULATION_WITH_MEMENTO_H___
#include <vector>
#include <algorithm>
using namespace std;

#include <szParticleF.h>
//#include <DisjointSet.h>

namespace Triangulation
{
	namespace _Internal
	{
		struct _triangle;
		struct _edge;

		bool isBoundary(const _edge* e);
		bool isBoundary(const _triangle* t);
		bool isActive(const _edge* e);
		bool isActive(const _triangle* e);

		struct _Memento
		{
			_edge* edge;
			_triangle* face;
		};

		struct _vertex
		{
			_vertex(const CParticleF& a): p(a)
			{
				bActive = true;
			}
			bool operator==(const _vertex& v) const
			{
				return p==v.p;
			}
			bool RemoveEdge(_edge* e)
			{
				//to remove it, e needs to be inactivated. nothing should change on vertex.
				return true;
			}
			bool AddEdge(_edge* e)
			{
				vector<_edge*>::iterator it = find(edges.begin(), edges.end(), e);
				if(it == edges.end())
				{
					edges.push_back(e);
					return true;
				}
				else
					return false;
			}
			CParticleF p;
			vector<_edge*> edges;
			bool bActive;
		};

		enum TriangulationEdgeType {Unknown, Inside, Hole, Boundary, Isolated};

		struct _edge
		{
			_edge(_vertex* a, _vertex* b)
			{
				vertices[0] = a;
				vertices[1] = b;
				faces[0] = NULL; faces[1]=NULL;
				length = Distance(vertices[0]->p, vertices[1]->p);
				type = Unknown;
				bActive = true;
			}
			bool Equal(const _vertex* a, const _vertex* b) const
			{
				return *vertices[0]==*a && *vertices[1]==*b ||
					*vertices[0]==*b && *vertices[1]==*a;
			}
			bool operator==(const _edge& e) const
			{
				return Equal(e.vertices[0], e.vertices[1]);
			}
			float Length() const
			{
				return length; //Distance(vertices[0]->p, vertices[1]->p);
			}
			bool operator <(const _edge& e) const
			{
				return Length() < e.Length();
			}
			bool AddFace(_triangle* t)
			{
				if(faces[0]==NULL) 
				{
					faces[0] = t;
					return true;
				}
				else if(faces[1]==NULL)
				{	
					faces[1] = t;
					return true;
				}
				else 
				{
					return false;
				}
			}
			bool RemoveFace(_triangle* t)
			{
				//to remove face, t has to be inactivated. nothing needs to be done on this edge.
				return true;
			}
			void SetLength(float len) {length = len;}

			_vertex* vertices[2];
			_triangle* faces[2];
			float length;
			TriangulationEdgeType type;
			bool bActive;
		};

		struct _edge_equal
		{
			_edge_equal(const _edge* e) : edge(e) { }
			bool operator()(_edge* e) const {return  *e == *edge;}  
		private:
			const _edge* edge;
		};		
		inline bool _edge_less(const _edge* e1, const _edge* e2) 
		{	
			return *e1 < *e2;
		}

		struct _triangle
		{
			_triangle(_edge* a, _edge* b, _edge* c, _vertex* u=0, _vertex* v=0, _vertex* w=0)
			{
				edges[0] = a; edges[1]=b; edges[2]=c;
				vertices[0]=u; vertices[1]=v; vertices[2]=w;
				bActive = true;
			}
			bool operator==(const _triangle& e) const
			{
				return *edges[0]==*e.edges[0] && *edges[1]==*e.edges[1] && *edges[2]==*e.edges[2] ||
					*edges[0]==*e.edges[1] && *edges[1]==*e.edges[2] && *edges[2]==*e.edges[0] ||
					*edges[0]==*e.edges[2] && *edges[1]==*e.edges[0] && *edges[2]==*e.edges[1];
			}
			_vertex* vertices[3];
			_edge* edges[3];
			bool bActive;
		};
		struct _triangle_equal{
			_triangle_equal(const _triangle* t) : tri(t) { }
			bool operator()(_triangle* t) const {return  *t == *tri;}  
		private:
			const _triangle* tri;
		};		

		/*	The result is <0 if the three points are in counter-clock wise.
		It is >0 if they are clock wise.
		It is =0 if they are colinear.
		*/

		struct _indexed_triangle
		{
			int index[3];
		};

	}

	int CounterClockWise(const CParticleF& p1, const CParticleF& p2, const CParticleF& p3);
	_Internal::_vertex* findOpposite(const _Internal::_triangle* tr, const _Internal::_edge* e);

	class Triangulator
	{
	public:
		Triangulator(const vector<CParticleF>& points);
		Triangulator(const vector<CParticleF>& points, const vector<_Internal::_indexed_triangle>& faces){
			_init(points, faces);
		}
		void _init(const vector<CParticleF>& pnts, const vector<_Internal::_indexed_triangle>& faces);
		~Triangulator()
		{
			for(int i=0; i<faces.size(); ++i)
			{
				delete faces[i];
			}
			for(int i=0; i<edges.size(); ++i)
			{
				delete edges[i];
			}
			for(int i=0; i<points.size(); ++i)
			{
				delete points[i];
			}
		}

		_Internal::_edge* GetNext() const;
		bool RemoveEdge(_Internal::_edge* e); //only a boundary edge can be removed.
		int SanityCheck() const;
		//private:
		void _DeleteEdge(_Internal::_edge* e); //utility to delete an edge and clean up
		vector<_Internal::_vertex*> points;
		vector<_Internal::_triangle*> faces; 
		vector<_Internal::_edge*> edges;  
		vector<_Internal::_Memento> states;

		//Memento
		int getState() 
		{ 
			return states.size();
		}
		bool setState(int state)
		{
			for(int i=states.size()-1; i>=state; --i)
			{
				states[i].edge->bActive = true;
				if(states[i].face)
				{
					states[i].face->bActive = true;
				}
			}
			states.erase(states.begin()+state, states.end());
			return true;
		}
		int numActiveEdges() const
		{
			int count = 0;
			for(int i=0; i<edges.size(); ++i)
			{
				if(edges[i]->bActive) count++;
			}
			return count;
		}
	};
}

#endif /* ___TRIANGULATION_WITH_MEMENTO_H___ */