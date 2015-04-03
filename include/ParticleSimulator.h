#pragma once;
#include <MovingParticle.h>
#include <mex.h>
#include <Snapshot.h>

class ParticleSimulator
{
public:
	ParticleSimulator()
	{
		time = 0.0f;
	}
	bool Prepare(vector<CParticleF>& points, vector<pair<int, int>>& edges, float delta0 = 0.01);
	bool LoadParticles(vector<float>& state, const int* dims);
	bool Restore(vector<Snapshot>& snapshots);
	mxArray* SaveParticles();
	mxArray* SaveDoneEvents();
	mxArray* SaveConvexity();
	bool Simulate(float endtime = 10.0f, float delta = 0.1f, bool bdebug = false);
	float getTime() const { return time; }

	static vector<MovingParticle*>	initializePolygon(vector<Edge<CParticleF>*>& edges); //temporary

	vector<Snapshot> snapshots; //id-time pair.
	vector<EventStruct> doneEvents;
	vector<Snapshot> closedRegions;
private:
	bool initializePolygon(vector<MovingParticle*>& particles);
	bool _Restore(Snapshot& snapshot);
	float time;
};

