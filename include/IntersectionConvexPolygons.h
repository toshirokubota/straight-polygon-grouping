#ifndef __INTERSECTION_CONVEX_POLYGONS_H___
#define __INTERSECTION_CONVEX_POLYGONS_H___
#include <vector>
#include <szParticleF.h>

namespace _IntersectConvexPolygon 
{
	/*
	find an intersection between two line segments defined by [p1 - p2] and [q1 - q2].
	It returns a parameter s such that the intersection point is (1-s)*p1 + s*p2.
	*/
	std::pair<float,float>
		intersect(const CParticleF& p1, const CParticleF& p2, const CParticleF& q1, const CParticleF& q2);

	/*
	Both left and right are convex polygons. 
	This function finds the intersection of the two convex polygons, which is also convex.
	*/
	std::vector<CParticleF> 
	intersectConvexHulls(const std::vector<CParticleF>& left,
						const std::vector<CParticleF>& right);
}

#endif /* __INTERSECTION_CONVEX_POLYGONS_H___ */