#pragma once

#include <szParticleF.h>
#include <vector>
#include <Graph.h>
#include <mex.h>
#include <map>
using namespace std;

namespace StraightAxis
{
	struct MovingParticle;
	struct DeformingPolygon;

	enum EventType {UnknownEvent, CollisionEvent, SplitEvent, MergeEvent, IntervalEvent};

	struct EventStruct
	{
		EventStruct(float t=std::numeric_limits<float>::infinity(), EventType type=UnknownEvent, MovingParticle* p=NULL, MovingParticle* q=NULL);
		bool operator <(EventStruct& ev)
		{
			return t < ev.t;
		}
		float t; 
		EventType type;
		MovingParticle* p;
		MovingParticle* q;
		MovingParticle* r;
	};

	enum MovingParticleType {Unknown, Initial, Regular, Merge, Collide, Split, Axis, Dummy};
	struct MovingParticle
	{
		MovingParticle(CParticleF& p=CParticleF(0,0,0), MovingParticleType t=Unknown, float tm=0.0f)
		{
			this->p = p;
			p0 = p;
			v[0] = 0; v[1] = 0;
			next = NULL;
			prev = NULL;
			id = _id++;
			type = t;
			created = tm;
			time = created;
			bStable = true;
			bActive = true;
			bDirty = false;
			event = EventStruct(std::numeric_limits<float>::infinity(), UnknownEvent, this);
		}
		~MovingParticle()
		{
			next = NULL;
			prev = NULL;
			anscestors.clear();
			descendents.clear();
		}
		bool setVelocity();
		CParticleF move0(float t)
		{
			//t = t - created;
			return CParticleF(p0.m_X+t*v[0], p0.m_Y+t*v[1]);
		}
		CParticleF move(float t)
		{
			//t = t - created;
			return CParticleF(p.m_X+t*v[0], p.m_Y+t*v[1]);
		}
		void update(float t)
		{
			//t = t - created;
			if(bActive)
			{
				p = move(t);
				time += t;
			}
		}
		void print(char* tab="")
		{
			printf("%s%d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f\n",	tab, id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1]);
			//printf("%d: [%3.3f,%3.3f]->[%3.3f,%3.3f] @(%3.3f,%3.3f)\n",	id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1]);
			//printf("%3.3f %3.3f %3.3f %3.3f\n",	p.m_X, p.m_Y, p0.m_X, p0.m_Y);
		}
		CParticleF p0;
		CParticleF p;
		float v[2];
		MovingParticle* next;
		MovingParticle* prev;
		int id;
		MovingParticleType type;
		float created; //time this is created.
		float time;
		bool bStable; //false if the velocity is not accurate.
		bool bActive; //true if it is still moving.
		bool bDirty; //true when its event may not be up-to-date.
		EventStruct event;
		vector<MovingParticle*> anscestors;
		vector<MovingParticle*> descendents;
		static int _id;
	};

	inline bool compareByEvent(MovingParticle* p, MovingParticle* q)
	{
		return p->event < q->event;
	}

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
		MovingParticle* makeParticle(CParticleF& p, MovingParticleType type, float tm)
		{
			MovingParticle* particle = new MovingParticle(p, type, tm);
			particles.push_back(particle);
			if (particle->id == 62)
			{
				particle->id += 0;
			}
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
		DeformingPolygon(float tm=0, bool ori = false, DeformingPolygon* pa = NULL)
		{
			created = tm;
			time = tm;
			id = _id++;
			parent = pa;
			died = 0;
		}
		DeformingPolygon(vector<CParticleF>& pnts, float tm, DeformingPolygon* pa = NULL);
		DeformingPolygon(vector<MovingParticle*>& pnts, float tm, DeformingPolygon* pa = NULL);
		vector<DeformingPolygon*> deform(float time);
		vector<DeformingPolygon*> deform(float time, float interval);
		void removeUnstables();
		void _setMap()
		{
			pmap.clear();
			for(int i=0; i<particles.size(); ++i)
			{
				pmap[particles[i]] = i;
			}
		}
		void takeSnapshot()
		{
			snapshot.clear();
			for(int i=0; i<particles.size(); ++i)
			{
				snapshot.push_back(particles[i]->p);
			}
		}
		MovingParticle* applyInterval(MovingParticle* p, float delta);
		MovingParticle* applyCollision(MovingParticle* p, float delta);
		pair<MovingParticle*, MovingParticle*> applySplit(MovingParticle* p, MovingParticle* q, float delta);
		int sanityCheck();
		void quickFinish();

		float created; //time this polygon is created.
		float died; //time this polygon becomes inactive.
		float time;
		int id;
		DeformingPolygon* parent;
		vector<DeformingPolygon*> children;
		vector<MovingParticle*> particles;
		map<MovingParticle*,int> pmap;
		vector<CParticleF> snapshot;
		bool orientation; //true for out-growing (including the initial polygon).
		EventStruct event; //the event that produced this polygon.
		static int _id;
	};

	void
	setRelation(MovingParticle* child, MovingParticle* parent);

	vector<MovingParticle*>
		traceTree(vector<Vertex<CParticleF>*>& tree);

	vector<MovingParticle*>
		offsetPolygon(vector<MovingParticle*>& points, float margin, float eps);

	vector<MovingParticle*>
		traceBackPolygon(DeformingPolygon* polygon);

	/*
	For a trivial case, we can use this one to expedite processing.
	*/
	//void quickFinish(vector<MovingParticle*>& particles, float time);

	mxArray* StoreParticles(const vector<vector<MovingParticle*>>& polygons);
}
