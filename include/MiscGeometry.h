#pragma once
#include <szParticleF.h>
#include <vector>

/*
A unit vector from o to x in 2D
*/
CParticleF
NormalizedDirection(CParticleF& x, CParticleF& o);

/*
A unit vector from o to x in 3D
*/
CParticleF
NormalizedDirection3D(CParticleF& x, CParticleF& o);

/*
return a unit vector that bisects the angle formed by three points: a-o-b.
*/
CParticleF
bisector(CParticleF& o, CParticleF& a, CParticleF& b);

/*
Find a parameter t such that ol+t*u intersects a plane cutting a point op and perpendicular to z.
*/
float
intersectPlaneAndLine(CParticleF op, CParticleF& z, CParticleF& ol, CParticleF u);

CParticleF
crossProduct(CParticleF& a, CParticleF& b);

float
dotProduct(CParticleF& a, CParticleF& b);

/*
Given three points in 3D, define a unit vector perpendicular to the plane going through the three points.
The vector is always upward (i.e. its Z-coordinate is always positive). But this should not really matter...
*/
CParticleF
perp2Plane(CParticleF& a, CParticleF& b, CParticleF& c);

/*
Find the centroid of a polygon, POLY.
*/
CParticleF
centroidPoint(std::vector<CParticleF>& poly);

/*
Check if a-b-c lie on a line.
*/
bool
coLinear(CParticleF& a, CParticleF& b, CParticleF& c, float precision = 1.0e-3);

/*
Check if lines a-b and c-d are parallel.
*/
bool
isParallel(CParticleF& a, CParticleF& b, CParticleF& c, CParticleF& d, float precision = 1.0e-3);

/*
Check if a moving side (p-q) and a moving particle (a with velocity (vx,vy)) are approaching to each other.
*/
bool
isApproaching(CParticleF& a, float vx, float vy, CParticleF& p, CParticleF& q);

/*
Compute a unit vector normal to the side (p->q). The direction is to the left when we are at p and looking toward q.
*/
CParticleF
perpDirection(CParticleF& p, CParticleF& q);
