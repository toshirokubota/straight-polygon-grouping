#ifndef ___ANGULARLY_ORDERED_GRAPH_H___
#define ___ANGULARLY_ORDERED_GRAPH_H___
#include <Graph.h>
#include <ContourEQW.h>
#include <FragmentInfo.h>
#include <szParticleF.h>
#include <vector>
using namespace std;

class AngularlyOrderedGraph
{
public:
	vector<Vertex<FragmentInfo>*> vertices;
	CParticleF click;
	bool orientation;
	int iter;
	Vertex<FragmentInfo>* source;
	Vertex<FragmentInfo>* dest;

	AngularlyOrderedGraph() {orientation = true; iter=0; source=NULL; dest=NULL;}
	void ReverseDirection() {orientation = orientation ? false: true;}
	void Initialize(vector<Vertex<FragmentInfo>*>& v, vector<ContourEQW>& contours, FragmentInfo& start, FragmentInfo& end);
	//void SetInitialDirection(FragmentInfo& start, FragmentInfo& dest, vector<ContourEQW>& contours);
	bool Update(vector<ContourEQW>& contours, vector<FragmentInfo>& excluded, float (*UPDATE)(float, float));
	vector<Vertex<FragmentInfo>*> findSolution();
};

#endif /* ___ANGULARLY_ORDERED_GRAPH_H___ */