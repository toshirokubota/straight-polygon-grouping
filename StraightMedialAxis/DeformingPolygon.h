#pragma once

#include <szParticleF.h>
#include <vector>
#include <Graph.h>
using namespace std;

namespace StraightAxis
{
	enum MovingParticleType {Unknown, Regular, Merge, Collide, Split, Axis, Dummy};
	struct MovingParticle
	{
		MovingParticle(CParticleF& p=CParticleF(0,0,0), MovingParticleType t=Unknown)
		{
			this->p = p;
			p0 = p;
			v[0] = 0; v[1] = 0;
			next = NULL;
			prev = NULL;
			id = _id++;
			type = t;
		}
		bool setVelocity();
		CParticleF move(float t)
		{
			return CParticleF(p.m_X+t*v[0], p.m_Y+t*v[1]);
		}
		CParticleF p0;
		CParticleF p;
		float v[2];
		MovingParticle* next;
		MovingParticle* prev;
		int id;
		MovingParticleType type;
		static int _id;
	};

	struct ParticleFactory
	{
		static ParticleFactory* getInstance()
		{
			if(_instance == NULL)
			{
				_instance = new ParticleFactory();
			}
			return _instance;
		}
		MovingParticle* makeParticle(CParticleF& p, MovingParticleType type)
		{
			MovingParticle* particle = new MovingParticle(p, type);
			particles.push_back(particle);
			return particle;
		}
		void clean()
		{
			for(int i=0; i<particles.size(); ++i)
			{
				delete particles[i];
			}
			particles.clear();
			MovingParticle::_id = 0;
		}
		~ParticleFactory()
		{
			clean();
		}
		vector<MovingParticle*> particles;
	private:
		static ParticleFactory* _instance;
		ParticleFactory()
		{
			MovingParticle::_id = 0;
		}
	};

	struct DeformingPolygon
	{
		DeformingPolygon(){}
		DeformingPolygon(vector<CParticleF>& pnts);
		void deform(float t)
		{
			for(int i=0; i<particles.size(); ++i)
			{
				particles[i]->p = particles[i]->move(t);
			}
		}
		pair<float,MovingParticle*> nextEdgeEvent();
		pair<float,MovingParticle*> nextSplitEvent();
		//bool fixPolygon(vector<DeformingPolygon>& additions, float eps = 0.1f);
		//DeformingPolygon split(int begin, int end);
		int sanityCheck();
		vector<MovingParticle*> particles;
	};

	vector<DeformingPolygon>
		fixPolygon(DeformingPolygon& poly, float eps, bool debugOn);

	vector<CParticleF>
		traceTree(vector<Vertex<CParticleF>*>& tree);

	vector<CParticleF>
		offsetPolygon(vector<CParticleF>& points, float margin, float eps);

	/*
	For a trivial case, we can use this one to expedite processing.
	*/
	void
		quickFinish(vector<MovingParticle*> particles);
}
