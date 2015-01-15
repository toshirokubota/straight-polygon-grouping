#pragma once
#include <vector>
using namespace std;
#include <szParticleF.h>
#include <MovingParticle.h>
#include <mex.h>

class Snapshot
{
public:
	Snapshot(float t=0.0f)
	{
		time = t;
	}
	Snapshot(float t, vector<MovingParticle*>& particles)
	{
		time = t;
		for (int i = 0; i < particles.size(); ++i)
		{
			add(particles[i]->getId());
		}
	}
	static vector<Snapshot> TakeSnapshot(float time);
	static Snapshot LoadSnapshot(const mxArray* ptr);
	static vector<Snapshot> LoadSnapshots(const mxArray* ptr);
	static mxArray* StoreSnapshot(Snapshot& snapshot);
	static mxArray* StoreSnapshots(vector<Snapshot>& snapshot);
	int size() const { return particles.size(); }
	float getTime() const { return time; }
	int get(int i) const { return particles[i]; }
	bool operator==(const Snapshot& shot) const
	{
		return shot.size() == this->size() && shot._pset == this->_pset;
	}
private:
	void add(int id)
	{
		particles.push_back(id);
		_pset.insert(id);
	}
	float time;
	vector<int> particles;
	set<int> _pset;
};
