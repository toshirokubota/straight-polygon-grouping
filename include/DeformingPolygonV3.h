#pragma once

#include <szParticleF.h>
#include <vector>
#include <Graph.h>
#include <mex.h>
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
		CParticleF p0;
		CParticleF p;
		float v[2];
		MovingParticle* next;
		MovingParticle* prev;
		int id;
		MovingParticleType type;
		float created; //time this is created.
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
			if(particle->id==1164)
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

	struct DeformingPolygon
	{
		DeformingPolygon(float tm, int pa_id = -1)
		{
			created = tm;
			id = _id++;
			parent_id = pa_id;
		}
		DeformingPolygon(vector<CParticleF>& pnts, float tm, int pa_id = -1);
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
		float created; //time this polygon is created.
		int id;
		int parent_id;
		vector<MovingParticle*> particles;
		static int _id;
	};

	void
	setRelation(MovingParticle* child, MovingParticle* parent);

	vector<DeformingPolygon>
		fixPolygon(DeformingPolygon& poly, float time, float eps, bool debugOn);

	vector<MovingParticle*>
		traceTree(vector<Vertex<CParticleF>*>& tree);

	vector<MovingParticle*>
		offsetPolygon(vector<MovingParticle*>& points, float margin, float eps);

	/*
	For a trivial case, we can use this one to expedite processing.
	*/
	void quickFinish(vector<MovingParticle*>& particles, float time);

	int hierarcyConsistencyCheck();

	mxArray* StoreParticles(const vector<vector<MovingParticle*>>& polygons);

	mxArray* StorePolygons(vector<DeformingPolygon>& polygons);

	bool LoadPolygons(vector<DeformingPolygon>& polygons, const mxArray *prhs);

	mxArray* StoreStraightAxes(vector<DeformingPolygon>& polygons);

	bool LoadStraightAxes(vector<DeformingPolygon>& polygons, const mxArray *prhs);
}
