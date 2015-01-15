#ifndef ___BIDIRECTIONALDP_UTILES_H___
#define ___BIDIRECTIONALDP_UTILES_H___

#include <ContourEQW.h>
#include <Graph.h>
//#include <GraphFactory.h>
#include <FragmentInfo.h>

void
PrintCycle(vector<Vertex<FragmentInfo>*>& cycle);

/*
Stitch a number of EQW contours together into one EQW contour.
The orientation of input contours are assumed to be arranged properly in
by splitContours().
*/
ContourEQW
stitchContour(vector<Vertex<FragmentInfo>*>& cycle, 
			  vector<ContourEQW>& contours);

/*
Collect a number of contours according to the order given by the cycle.
It keeps multiple contours, not like stitchContour in which the contours
are combined into one.
*/
vector<ContourEQW>
collectContours(vector<Vertex<FragmentInfo>*>& cycle, 
			  vector<ContourEQW>& contours);

/*
Split contours so that each length is at most the second argument (length).
*/
vector<ContourEQW> splitContours(vector<ContourEQW>& contours, 
								 int length);

/*
split contours when the saliency value cross the threshold boundary
*/
vector<ContourEQW> 
splitContoursBySaliency(vector<ContourEQW>& contours, 
						float thres = 0.5);

/*
In addition to the length requirement, this version also splits a contour at the point closest to the start position.
*/
vector<ContourEQW> 
splitContoursByLengthAndAtSource(vector<ContourEQW>& contours0, 
								int maxLen, 
								CParticleF& start);

/*
Convert a path given by the first argument (vertices) and the second (contours) into a matrix and
return the result in the third argument.
It also returns the number of columns as the result in the fourth argument (ncols).
*/
void
Path2Mat(vector<Vertex<FragmentInfo>*>& vertices,
		  vector<ContourEQW>& contours,
		  vector<float>& res, 
		  int& ncols);

/*
Generate a set of vertices from a set of contours.
The number of vertices is twice the number of contours (one for the first and the other for the last
point in a contour).
*/
vector<Vertex<FragmentInfo>*>
generateVertices(vector<ContourEQW>& contours);

/*
Order vertices (1st arg) based on the visual angle seen from a click point (3rd arg).
The 2nd argument carries contour information.
*/
vector<Vertex<FragmentInfo>*>
orderVertices(vector<Vertex<FragmentInfo>*>& vertices,
			  vector<ContourEQW>& contours, 
			  CParticleF& click);

/*
Circularly shift vertices so that the source becomes the first one. 
Add a copy of the source at the end as the destination vertex.
*/
vector<Vertex<FragmentInfo>*>
arrangeVertices(vector<Vertex<FragmentInfo>*>& vertices,
				Vertex<FragmentInfo>* source);

/*
Given a position (x, y) as the 3rd and 4th arguments, find a contour point
that is the closest, and returns the vertex that is closest to the point.
*/
Vertex<FragmentInfo>* 
findSource(vector<Vertex<FragmentInfo>*>& vertices, 
		   vector<ContourEQW>& contours, 
		   float x, float y);

/*
Perform update on an angularly ordered graph (given as the 1st argument). The 2nd argument
carries the necessary contour information.
The third argument dictates the angular direction (true for contour clock wise, and false for
clock wise). The 4th argument restrict the difference in the view angle for two vertices to be linked.
Typically, this is set to pi/2.
Currently, the finalized time is set to the iteration number of the update. Thus, the 5th argument
is necessary.
*/
bool 
UpdateDistance(vector<Vertex<FragmentInfo>*>& vertices, vector<ContourEQW>& contours, 
			   bool orientation,
			   float angleThres,
			   int iter);

/*
A slight modification of UpdateDistance.
This version allows non-zero weights on solid fragments. The feature is used to give preference to fragments that are selected by the user.
It is used in BidirectionalDPMPv3.cpp.
*/
bool 
UpdateDistanceWithPreference(vector<Vertex<FragmentInfo>*>& vertices, vector<ContourEQW>& contours, 
							 bool orientation,
							 int iter);

/*
This is a slight modification to UpdateDistance. This version can exclude some nodes from being used
in the shortest path.
The fourth argument (excluded) supplies a set of vertices that need to be excluded from the update.
*/
bool 
UpdateDistanceWithExclusions(vector<Vertex<FragmentInfo>*>& vertices, vector<ContourEQW>& contours, 
			   bool orientation,
			   int iter,
			   vector<FragmentInfo>& excluded);

/*
Compute the total edge weights of the path.
*/
float
getPathCost(vector<Vertex<FragmentInfo>*>& path,
			vector<ContourEQW>& contours);

/*
Given an angularly ordered graph, trace a path from the 2nd argument (dest).
*/
vector<Vertex<FragmentInfo>*> 
findPath(vector<Vertex<FragmentInfo>*>& vertices, 
		  Vertex<FragmentInfo>* dest);

/*
Given an angularly ordered tree (vertices), find leaf vertices.
*/
vector<Vertex<FragmentInfo>*>
findEndPoints(vector<Vertex<FragmentInfo>*>& vertices);

/*
Given an angularly ordered tree (vertices), construct a set of contours where
each contour describe a branch of the tree.
*/
vector<ContourEQW> 
extractPathTree(vector<Vertex<FragmentInfo>*>& vertices, 
				vector<ContourEQW>& contours);


#endif /* ___BIDIRECTIONALDP_UTILES_H___ */