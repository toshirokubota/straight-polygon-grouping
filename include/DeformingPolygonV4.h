#pragma once

#include <szParticleF.h>
#include <vector>
#include <Graph.h>
#include <mex.h>
#include <map>
using namespace std;

namespace StraightAxis
{
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
		}
		~MovingParticle()
		{
			next = NULL;
			prev = NULL;
			anscestors.clear();
			descendents.clear();
		}
		bool setVelocity();
		CParticleF move(float t)
		{
			//t = t - created;
			return CParticleF(p.m_X+t*v[0], p.m_Y+t*v[1]);
		}
		void update(float t)
		{
			//t = t - created;
			p = move(t);
			time += t;
		}
		void print()
		{
			printf("%d: [%3.3f,%3.3f]->[%3.3f,%3.3f] @(%3.3f,%3.3f)\n",	id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1]);
			//printf("%3.3f %3.3f\n",	p.m_X, p.m_Y);
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
		vector<MovingParticle*> anscestors;
		vector<MovingParticle*> descendents;
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
		MovingParticle* makeParticle(CParticleF& p, MovingParticleType type, float tm)
		{
			MovingParticle* particle = new MovingParticle(p, type, tm);
			particles.push_back(particle);
			if(particle->id==1056)
				p.m_Life += 0;

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

	enum EventType {UnknownEvent, CollisionEvent, SplitEvent, MergeEvent};

	struct EventStruct
	{
		EventStruct(float t=0.0f, EventType type=UnknownEvent, MovingParticle* p=NULL, MovingParticle* q=NULL)
		{
			this->t = t;
			this->type = type;
			this->p = p;
			this->q = q;
			if(q != NULL && type == SplitEvent) //remember the split event, since this neighbor relation will be destroyed by the event
			{
				r = q->next;
			}
			else
			{
				r = NULL;
			}
		}
		bool operator <(EventStruct& ev)
		{
			return t < ev.t;
		}
		void print()
		{
			if(type==CollisionEvent)
			{
				printf("(%3.3f) C: [%3.3f,%3.3f]%d -> [%3.3f, %3.3f]%d\n", 
					t, p->p.m_X, p->p.m_Y, p->id, q->p.m_X, q->p.m_Y, q->id);
			}
			else if(type==SplitEvent)
			{
				printf("(%3.3f) S: [%3.3f,%3.3f]%d -> [%3.3f, %3.3f]%d - [%3.3f, %3.3f]%d\n", 
					t, p->p.m_X, p->p.m_Y, p->id, q->p.m_X, q->p.m_Y, q->id, r->p.m_X, r->p.m_Y, r->id);
			}
		}
		float t; 
		EventType type;
		MovingParticle* p;
		MovingParticle* q;
		MovingParticle* r;
	};

	struct DeformingPolygon
	{
		DeformingPolygon(float tm=0, bool ori = false, int pa_id = -1)
		{
			created = tm;
			time = tm;
			id = _id++;
			parent_id = pa_id;
		}
		DeformingPolygon(vector<CParticleF>& pnts, float tm, int pa_id = -1);
		DeformingPolygon(vector<MovingParticle*>& pnts, float tm, int pa_id = -1);
		vector<DeformingPolygon> deform(float delta);
		void _setMap()
		{
			pmap.clear();
			for(int i=0; i<particles.size(); ++i)
			{
				pmap[particles[i]] = i;
			}
		}
		MovingParticle* applyCollision(MovingParticle* p, float delta);
		pair<MovingParticle*,MovingParticle*> applySplit(MovingParticle* p, MovingParticle* q, float delta);
		int sanityCheck();
		float created; //time this polygon is created.
		float time;
		int id;
		int parent_id;
		vector<MovingParticle*> particles;
		map<MovingParticle*,int> pmap;
		bool orientation; //true for out-growing (including the initial polygon).
		static int _id;
	};

	void
	setRelation(MovingParticle* child, MovingParticle* parent);

	vector<DeformingPolygon>
		fixPolygon(DeformingPolygon& poly, vector<EventStruct>& events, float time, bool debugOn);

	vector<MovingParticle*>
		traceTree(vector<Vertex<CParticleF>*>& tree);

	vector<MovingParticle*>
		offsetPolygon(vector<MovingParticle*>& points, float margin, float eps);

	/*
	For a trivial case, we can use this one to expedite processing.
	*/
	void quickFinish(vector<MovingParticle*>& particles, float time);

	mxArray* StoreParticles(const vector<vector<MovingParticle*>>& polygons);
}
