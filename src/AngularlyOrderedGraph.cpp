#include <AngularlyOrderedGraph.h>
#include <BidirectionalDPUtils.h>
#include <szMiscOperations.h>
#include <szMexUtility.h>

bool GetOrientation(FragmentInfo& start, 
					FragmentInfo& dest,
					CParticleF& click,
					vector<ContourEQW>& contours)
{
	float x1 = contours[start.contourID].X[start.pointID];
	float y1 = contours[start.contourID].Y[start.pointID];
	float x2 = contours[dest.contourID].X[dest.pointID];
	float y2 = contours[dest.contourID].Y[dest.pointID];
	float d1 = GetVisualDirection(x1, y1, click.m_X, click.m_Y);
	float d2 = GetVisualDirection(x2, y2, click.m_X, click.m_Y);
	if(d2 < d1)
	{
		d2 += 2*PI;
	}
	if(d2-d1 < PI) return true;
	else return false;
}

void 
AngularlyOrderedGraph::Initialize(vector<Vertex<FragmentInfo>*>& v, 
								  vector<ContourEQW>& contours, 
								  FragmentInfo& start,
								  FragmentInfo& end)
{
	for(int i=0; i<vertices.size(); ++i)
	{
		delete vertices[i];
	}
	vertices.clear();

	for(int i=0; i<v.size(); ++i)
	{
		Vertex<FragmentInfo>* u = new Vertex<FragmentInfo>(v[i]->key);
		*u = *(v[i]);
		vertices.push_back(u);
	}

	vertices = orderVertices(vertices, contours, click);
	source = findVertex(vertices, start);
	assert(source);
	dest = findVertex(vertices, end);
	assert(dest);
	if(GetOrientation(start, end, click, contours)==false)
	{
		swap(source, dest);
	}

	vertices = arrangeVertices(vertices, source);

	for(int i=0; i<vertices.size(); ++i)
	{
		vertices[i]->Reset();
	}
	source->d = 0;
	iter = 1;
}

bool
AngularlyOrderedGraph::Update(vector<ContourEQW>& contours,
							  vector<FragmentInfo>& excluded,
							  float (*UPDATE)(float, float))
{
	bool b = UpdateDistanceWithExclusions(vertices, contours, orientation, iter, excluded, UPDATE);
	orientation = orientation ? false: true;
	iter++;
	return b;
}

vector<Vertex<FragmentInfo>*> 
AngularlyOrderedGraph::findSolution()
{
	Vertex<FragmentInfo>* v = findVertex(vertices, dest->key);
	assert(v);
	return findPath(vertices, v);
}


