#include <TriangulationHelper.h>
#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>

std::vector<CParticleF> 
vector2particle(const std::vector<float>& P,
				const int* dims)
{
	int nump = dims[0];
	vector<CParticleF> points(nump);
	for(int i=0; i<nump; ++i)
	{
		points[i].m_X = GetData2(P, i, 0, dims[0], dims[1], (float)0);
		points[i].m_Y = GetData2(P, i, 1, dims[0], dims[1], (float)0);
		/*if(dims[1] >= 3)
		{
			points[i].m_Z = GetData2(P, i, 2, dims[0], dims[1], (float)0);
		}
		if(dims[1] >= 4)
		{
			points[i].m_Life = GetData2(P, i, 3, dims[0], dims[1], (float)0);
		}*/
	}
	return points;
}

std::vector<Triangulation::_Internal::_indexed_triangle> 
indices2structs(const std::vector<CParticleF>& P, 
				const std::vector<int>& T, 
				const int* dimsT)
{
	int n = dimsT[0];
	vector<Triangulation::_Internal::_indexed_triangle> tr(n);
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			int k = GetData2(T, i, j, dimsT[0], dimsT[1], -1);
			if(k<1) abort();
			tr[i].index[j] = k - 1; //change from 1 to 0 index.
		}
	}
	return tr;
}

/*
get indices in [m, n] at increment of 'skip'. 
*/
vector<int>
getIndices(int m, int n, int skip)
{
	vector<int> indices;
	if(n==m) indices.push_back(m);
	else if(n > m)
	{
		indices.push_back(m);
		int k = floor((float)(n - m - 1) / (float)skip);
		if(k > 0)
		{
			float inc = (float)(n - m) / (k + 1);
			float z = m + inc;
			while(floor(z) < n)
			{
				indices.push_back(floor(z));
				z += inc;
			}
		}
		indices.push_back(n);
	}
	return indices;
}

/*
Select every 'interval' points
*/
vector<vector<FragmentInfo>>
TriangulationPointsFixedInterval(const vector<ContourEQW>& contours, int interval)
{
	vector<vector<FragmentInfo>> placements(contours.size()); //contour points selected for triangulation
	//first select all contour end-points.
	for(int i=0; i<contours.size(); ++i)
	{
		vector<int> range = getIndices(0, contours[i].size()-1, interval);
		for(int j=0; j<range.size(); ++j)
		{
			placements[i].push_back(FragmentInfo(i, range[j]));
		}
	}

	return placements;
}

/*
Sample uniformly so that the number of points is NUMP.
*/
vector<vector<FragmentInfo>>
TriangulationPointsSampleSize(const vector<ContourEQW>& contours, int nump)
{
	vector<vector<FragmentInfo>> placements(contours.size()); //contour points selected for triangulation
	//first select all contour end-points.
	for(int i=0; i<contours.size(); ++i)
	{
		int interval = contours[i].size() / nump;
		for(int j=0; j<contours[i].size(); j+=interval)
		{
			placements[i].push_back(FragmentInfo(i, j));
		}
	}

	return placements;
}

float
flatness(Triangulation::_Internal::_triangle* tr, Triangulation::_Internal::_edge* e)
{
	float minlen = std::numeric_limits<float>::infinity();
	for(int i=0; i<3; ++i)
	{
		if(tr->edges[i] != e)
		{
			minlen = Min(minlen, tr->edges[i]->Length());
		}
	}
	return e->Length() - minlen;
}

float EffectiveLength(Triangulation::_Internal::_edge* e)
{
	return e->Length();

	float len = 0;
	for(int i=0; i<2; ++i)
	{
		if(e->faces[i] != NULL)
		{
			float val = flatness(e->faces[i], e);
			if(val > len)
			{
				len = val;
			}
		}
	}
	return len;
}

vector<CParticleF>
traceShape(const vector<Triangulation::_Internal::_triangle*> faces)
{
	vector<CParticleF> shape;
	//pick the initial points
	Triangulation::_Internal::_vertex* cur = NULL;
	Triangulation::_Internal::_edge* ed = NULL;
	for(int i=0; i<faces.size() && cur==NULL; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			if(faces[i]->edges[j]->type == Triangulation::_Internal::Boundary)
			{
				ed = faces[i]->edges[j];
				cur = faces[i]->edges[j]->vertices[0];
				shape.push_back(cur->p);
				break;
			}
		}
	}

	while(cur)
	{
		Triangulation::_Internal::_vertex* next = ed->vertices[0]==cur ? ed->vertices[1]: ed->vertices[0];
		if(find(shape.begin(), shape.end(), next->p) != shape.end()) break;
		cur = next;
		shape.push_back(cur->p);
		for(int i=0; i<cur->edges.size(); ++i)
		{
			if(cur->edges[i] != ed && cur->edges[i]->type == Triangulation::_Internal::Boundary)
			{
				ed = cur->edges[i];
				break;
			}
		}
	}
	return shape;
}

bool
isExposed(Triangulation::_Internal::_vertex* v)
{
	for(int j=0; j<v->edges.size(); ++j)
	{
		if(v->edges[j]->type == Triangulation::_Internal::Boundary)
		{
			return true;
		}
	}
	return false;
}

Triangulation::_Internal::_vertex*
oppositeSideVertex(Triangulation::_Internal::_edge* e, Triangulation::_Internal::_triangle* tr)
{
	assert(e==tr->edges[0] || e==tr->edges[1] || e==tr->edges[2]); //make sure e is a part of tr.
	for(int i=0; i<3; ++i)
	{
		if(tr->vertices[i] != e->vertices[0] && tr->vertices[i] != e->vertices[1])
		{
			return tr->vertices[i];
		}
	}
	return NULL;
}

/*
find the vertex opposite to a BOUNDARY edge, e.
*/
Triangulation::_Internal::_vertex* 
oppositeSideVertexOfBoundaryEdge(Triangulation::_Internal::_edge* e)
{
	Triangulation::_Internal::_triangle* face = NULL;
	if(e->faces[0]!=NULL && e->faces[1]==NULL)
	{
		face = e->faces[0];
	}
	else if(e->faces[0]==NULL && e->faces[1]!=NULL)
	{
		face = e->faces[1];
	}
	else
	{
		return NULL; //this should not happen...
	}
	return oppositeSideVertex(e, face);
}

/*
Find a veretx that is commont to the given two edges in a triangulated graph.
*/
Triangulation::_Internal::_vertex*
commonVertex(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2)
{
	Triangulation::_Internal::_vertex* u = NULL;
	for(int i=0; i<2 && u==NULL; ++i)
	{
		for(int j=0; j<2; ++j)
		{
			if(e1->vertices[i]==e2->vertices[j])
			{
				u = e1->vertices[i];
				break;
			}
		}
	}
	return u;
}

/*
Given a pair of vertices, fined an edge (if any) that connects the two vertices.
If no edge is present between them, return NULL.
*/
Triangulation::_Internal::_edge*
findEdge(Triangulation::_Internal::_vertex* u, Triangulation::_Internal::_vertex* v)
{
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* ed = u->edges[i];
		if((ed->vertices[0]==u && ed->vertices[1]==v) || (ed->vertices[0]==v && ed->vertices[1]==u))
		{
			return ed;
		}
	}
	return NULL;
}

/*
Find the edge that lies between the given two triangles.
*/
Triangulation::_Internal::_edge*
commonEdge(Triangulation::_Internal::_triangle* t, Triangulation::_Internal::_triangle* s)
{
	for(int i=0; i<3; ++i)
	{
		if(t->edges[i]->faces[0]==s || t->edges[i]->faces[1]==s)
		{
			return t->edges[i];
		}
	}
	return NULL;
}

/*
Find the anglee (in radian) between two edges, which have to be incident to each other.
*/
float 
angleBetween(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2)
{
	Triangulation::_Internal::_vertex* u = commonVertex(e1, e2);
	Triangulation::_Internal::_vertex* v = e1->vertices[0] == u ? e1->vertices[1]: e1->vertices[0];
	Triangulation::_Internal::_vertex* w = e2->vertices[0] == u ? e2->vertices[1]: e2->vertices[0];
	return GetVisualAngle(v->p.m_X, v->p.m_Y, w->p.m_X, w->p.m_Y, u->p.m_X, u->p.m_Y);
}

/*
Remove edges without faces (dangling edges) from the triangulation.
*/
void
RemoveDanglingEdges(Triangulation::Triangulator& trmap)
{
	for(int i=trmap.edges.size()-1; i>=0; --i)
	{
		Triangulation::_Internal::_edge* e = trmap.edges[i];
		if(e->faces[0] == NULL && e->faces[1] == NULL)
		{
			trmap.RemoveIsolatedEdge(e);
		}
	}
}
