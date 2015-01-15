#include <ConvexityMeasure.h>
//#include <Triangulation.h>
#include <TriangulationHelper.h>
#include <szMiscOperations.h>

float
ConvexityMeasure(const std::vector<CParticleF>& points, float maxGap)
{
	Triangulation::Triangulator trmap(points);
	while(true)
	{
		vector<Triangulation::_Internal::_edge*> toRemove;
		for(int i=0; i<trmap.edges.size(); ++i)
		{
			if(trmap.edges[i]->type == Triangulation::_Internal::Boundary) 
			{
				if(trmap.edges[i]->Length() > maxGap)
				{
					toRemove.push_back(trmap.edges[i]);
				}
			}
		}
		if(toRemove.empty()) break;
		for(int j=0; j<toRemove.size(); ++j)
		{
			trmap.RemoveAnyEdge(toRemove[j]);
		}
		RemoveDanglingEdges(trmap);
	}
	float measure = ConvexityMeasure(trmap);
	return measure;
}

float
ConvexityMeasure(Triangulation::Triangulator& trmap)
{
	float sum = 0.0f;
	for(int i=0; i<trmap.deleted_faces.size(); ++i)
	{
		sum += areaTriangle(trmap.deleted_faces[i]->vertices[0]->p, trmap.deleted_faces[i]->vertices[1]->p, trmap.deleted_faces[i]->vertices[2]->p);
	}

	float sum2 = 0.0f;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		sum2 += areaTriangle(trmap.faces[i]->vertices[0]->p, trmap.faces[i]->vertices[1]->p, trmap.faces[i]->vertices[2]->p);
	}
	
	return sum2 / (sum + sum2);
}
