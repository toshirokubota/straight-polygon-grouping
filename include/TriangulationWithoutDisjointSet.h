#ifndef ___TRIANGULATION_H___
#define ___TRIANGULATION_H___
#include <vector>
#include <algorithm>
using namespace std;

#include <szParticleF.h>
#include <DisjointSet.h>

namespace Triangulation
{
	namespace _Internal
	{
		struct _triangle;
		struct _edge;

		struct _vertex
		{
			_vertex(const CParticleF& a): p(a)
			{
			}
			bool operator==(const _vertex& v) const
			{
				return p==v.p;
			}
			bool RemoveEdge(_edge* e)
			{
				vector<_edge*>::iterator it = find(edges.begin(), edges.end(), e);
				if(it != edges.end())
				{
					edges.erase(it);
					return true;
				}
				else
					return false;
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
				if(faces[1] == t) 
				{
					faces[1] = NULL;
					return true;
				}
				else if(faces[0] == t) 
				{
					//faces[0] = faces[1];
					faces[0] = NULL;
					return true;
				}
				else return false;
			}
			void SetLength(float len) {length = len;}

			_vertex* vertices[2];
			_triangle* faces[2];
			float length;
			TriangulationEdgeType type;
		};

		struct _edge_equal{
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
			}
			bool operator==(const _triangle& e) const
			{
				return *edges[0]==*e.edges[0] && *edges[1]==*e.edges[1] && *edges[2]==*e.edges[2] ||
					*edges[0]==*e.edges[1] && *edges[1]==*e.edges[2] && *edges[2]==*e.edges[0] ||
					*edges[0]==*e.edges[2] && *edges[1]==*e.edges[0] && *edges[2]==*e.edges[1];
			}
			_vertex* vertices[3];
			_edge* edges[3];
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
		Triangulator(const vector<CParticleF>& points, const vector<_Internal::_indexed_triangle>& faces);

		~Triangulator()
		{
			for(int i=0; i<faces.size(); ++i)
			{
				delete faces[i];
			}
			for(int i=0; i<deleted_faces.size(); ++i)
			{
				delete deleted_faces[i];
			}
			for(int i=0; i<edges.size(); ++i)
			{
				delete edges[i];
			}
			for(int i=0; i<deleted_edges.size(); ++i)
			{
				delete deleted_edges[i];
			}
			for(int i=0; i<points.size(); ++i)
			{
				delete points[i];
			}
		}
		_Internal::_edge* GetNext() const;
		bool RemoveSafeEdge(_Internal::_edge* e); //this used to be RemoveEdge. It only removes a safe edge
		bool RemoveAnyEdge(_Internal::_edge* e); //this remove any edge, safe or not.
		bool RemoveIsolatedEdge(_Internal::_edge* e); //this remove a dangling edge (i.e. an edge without triangles attached to it).
		int SanityCheck() const;
		//private:
		void _DeleteEdge(_Internal::_edge* e); //utility to delete an edge and clean up
		vector<_Internal::_vertex*> points;
		vector<_Internal::_triangle*> faces; 
		vector<_Internal::_edge*> edges;  
		vector<_Internal::_triangle*> deleted_faces; 
		vector<_Internal::_edge*> deleted_edges;  
	};
}

#endif /* ___TRIANGULATION_H___ */