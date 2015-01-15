#ifndef ___STRAIGHT_OFFSET_POLYGONS_H___
#define ___STRAIGHT_OFFSET_POLYGONS_H___
#include <vector>
#include <set>
#include <limits>
using namespace std;
#include <szParticleF.h>
#include <MiscGeometry.h>
#include <MovingParticle.h>
#include <Graph.h>

class StraightOffsetPolygon
{
public:
	StraightOffsetPolygon(const vector<MovingParticle*>& vertices, float time);
	pair<StraightOffsetPolygon*, StraightOffsetPolygon*> applyNextEvent();
	vector<CParticleF> snapShot() const;
	vector<MovingParticle*> traceBack() const;
	MovingParticle* get(int i) { return particles[i]; }
	int size() const { return particles.size(); }
	bool exist(const MovingParticle* p) const
	{
		return pset.find(p) != pset.end();
	}
private:
	vector<MovingParticle*> particles;
	set<const MovingParticle*> pset;
	float created;
};

/*
Utility functions
*/
/*
Trace a forrest and create a polygon for each tree.
*/
vector<StraightOffsetPolygon*>
traceForrest(vector<Vertex<CParticleF>*>& forrest);

#endif /* ___STRAIGHT_OFFSET_POLYGONS_H___*/