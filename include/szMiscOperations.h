#ifndef _SZ_MISC_OPERATIONS_H_
#define _SZ_MISC_OPERATIONS_H_
#include <vector>
using namespace std;
#include <szParticleF.h>
const double PI = 3.1415926536;

int GetOffsetX6(int n);

int GetOffsetY6(int n);

int GetOffsetZ6(int n);

int GetOffsetX10(int n);

int GetOffsetY10(int n);

int GetOffsetZ10(int n);

int GetOffsetX26(int n);

int GetOffsetY26(int n);

int GetOffsetZ26(int n);

/*
Get the median of the three numbers.
*/
int 
Middle(int a, int b, int c);

//float
//Distance(const CParticleF& p1, const CParticleF& p2);

float
Distance(float x1, float y1, float x2, float y2);

/*
Check if a trace of three points (x0,y0) -> (x1,y1) -> (x2,y2) is in Clockwise.
It returns 1 if clockwise, -1 if counter clockwise, and 0 if undetermined.
*/
int
ClockWise(float x0, float y0, float x1, float y1, float x2, float y2);

/*
Check if a trace of points is in Clockwise.
It returns 1 if clockwise, -1 if counter clockwise, and 0 if undetermined.
*/
int
ClockWise(const vector<CParticleF>& points);

/*
get the direction of (x, y) as seen from (x0, y0).
*/
float
GetVisualDirection(float x, float y, 
				   float x0, float y0);

/*
Get the angle of view formed by (x1, y1) and (x2, y2) as seen from (x0, y0).
*/
float
GetVisualAngle(float x1, float y1, float x2, float y2, 
			   float x0, float y0);

/*
It considers the orientation of (x1, y1) and (x2, y2) w.r.t. (x0,y0).
It returns negative if (x2, y2) is not clockwise orientation to (x1, y1).
It returns positive otherwise.
*/
float
GetVisualAngle2(float x1, float y1, float x2, float y2, 
			   float x0, float y0);

CParticleF
representativePosition(const vector<CParticleF>& cluster);

float
areaTriangle(CParticleF& a, CParticleF& b, CParticleF& c);

bool
inside(const CParticleF& p, const vector<CParticleF>& pnts);

/*
compute the closest point from x on the line connecting p and q.
p and q have to be distinct.
*/
CParticleF Closest2Line(const CParticleF& p, const CParticleF& q, const CParticleF& x);

/*
compute the distance from x to the line connecting a and b.
*/
float Distance2Line(const CParticleF& a, const CParticleF& b, const CParticleF& x);

float Distance2LineSegment(const CParticleF& p, const CParticleF& q, const CParticleF& x);

float
polygonArea(vector<CParticleF>& vp);

#endif /* _SZ_MISC_OPERATIONS_H_ */