#include <BidirectionalDPUtils.h>
#include <GraphFactory.h>

/*
Generate a set of vertices from a set of contours.
The number of vertices is twice the number of contours (one for the first and the other for the last
point in a contour).
*/
vector<Vertex<FragmentInfo>*>
generateVertices(vector<ContourEQW>& contours)
{
	GraphFactory<FragmentInfo>* factory = GraphFactory<FragmentInfo>::GetInstance();
	vector<Vertex<FragmentInfo>*> vertices;
	for(int i=0; i<contours.size(); ++i)
	{
		Vertex<FragmentInfo>* u = factory->makeVertex(FragmentInfo(i, 0));
		Vertex<FragmentInfo>* v = factory->makeVertex(FragmentInfo(i, contours[i].size()-1));
		vertices.push_back(u);
		vertices.push_back(v);
	}
	return vertices;
}


