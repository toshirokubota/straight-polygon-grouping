#ifndef ___TRIANGULATE_BOURKE_H___
#define ___TRIANGULATE_BOURKE_H___
#include <vector>
#include <szParticleF.h>

namespace TriangulateBourke
{
	namespace _Internal
	{
		typedef struct 
		{
			int p1,p2,p3;
		} ITRIANGLE;
		typedef struct 
		{
			int p1,p2;
		} IEDGE;
		typedef struct 
		{
			double x,y,z;
		} XYZ;
	}

	struct Triplet
	{
		int index[3];
	};

	std::vector<Triplet> DelaunayTriangulation(const std::vector<CParticleF>& points);

}

#endif /*  ___TRIANGULATE_BOURKE_H___ */