#pragma once
#include <vector>
#include <set>
using namespace std;

#include <szParticleF.h>
#include <EventStruct.h>
#include <MiscGeometry.h>

enum MovingParticleType { Unknown, Initial, Regular, Merge, Collide, Split, Axis, Dummy };
struct ParticleFactory;

class MovingParticle
{
public:
	bool setVelocity(float vx, float vy)
	{
		if (bInitialized == false)
		{
			//velocity can be set only once.
			v[0] = vx;
			v[1] = vy;
			bInitialized = true;
			return true;
		}
		else
		{
			return false;
		}
	}
	void update(float delta)
	{
		//if (bActive)
		{
			p = move(delta);
			time += delta;
		}
	}
	void print(char* tab = "") const
	{
		printf("%s%d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f\n", tab, id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1]);
	}

	CParticleF getP0() const { return p0; }
	CParticleF getP() const { return p; }
	void getVelocity(float& vx, float& vy) { vx = v[0]; vy = v[1]; }
	MovingParticle* getNext() const { return next; }
	MovingParticle* getPrev() const { return prev; }
	int getId() const { return id; }
	EventStruct getEvent() const { return event; }
	float getTime() const { return time; }
	MovingParticleType getType() const { return type; }
	//bool isActive() const { return bActive; }
	bool isInitialized() const { return bInitialized; }

	static void setNeighbors(MovingParticle* p, MovingParticle* prev, MovingParticle* next)
	{
		p->prev = prev;
		p->next = next;
		p->prev->next = p;
		p->next->prev = p;
	}
	static vector<MovingParticle*> vectorize(MovingParticle* p);
	bool updateEvent();
	pair<MovingParticle*, MovingParticle*> applyEvent();

	friend ParticleFactory; //so that the factory can access _id.
private:
	MovingParticle(CParticleF& p = CParticleF(0, 0, 0), MovingParticleType t = Unknown, float tm = 0.0f)
	{
		this->p = p;
		p0 = p;
		v[0] = 0; v[1] = std::numeric_limits<float>::quiet_NaN();
		next = NULL;
		prev = NULL;
		anscestors[0] = anscestors[1] = NULL;
		descendents[0] = descendents[1] = NULL;
		id = _id++;
		type = t;
		created = tm;
		time = created;
		//bActive = true;
		bInitialized = false;
		event = EventStruct(std::numeric_limits<float>::infinity(), UnknownEvent, this);
	}
	CParticleF move(float t) const
	{
		return CParticleF(p.m_X + t*v[0], p.m_Y + t*v[1]);
	}
	static float intersectPolygonAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r);
	EventStruct findNextEdgeEvent() const;
	EventStruct findNextSplitEvent() const;

	CParticleF p0;
	CParticleF p;
	MovingParticle* next;
	MovingParticle* prev;
	MovingParticle* anscestors[2];
	MovingParticle* descendents[2];
	int id;
	MovingParticleType type;
	float created; //time this is created.
	float time;
	bool bInitialized; //false until its velocity is set.
	//bool bActive; //true if it is still moving.
	EventStruct event;
	float v[2];
	static int _id;
};

/*struct ParticleQueue
{
	static ParticleQueue* getInstance()
	{
		if (_instance == NULL)
		{
			_instance = new ParticleQueue();
		}
		return _instance;
	}
	void add(MovingParticle* p)
	{
		pset.insert(p);
	}
	void  remove(MovingParticle* p)
	{
		pset.erase(p);
	}
	set<MovingParticle*>::const_iterator iterator()
	{
		return pset.begin();
	}
	set<MovingParticle*>::const_iterator end()
	{
		return pset.end();
	}

private:
	static ParticleQueue* _instance;
	ParticleQueue()
	{
	}
	set<MovingParticle*> pset;
};*/

struct ParticleFactory
{
	static ParticleFactory* getInstance()
	{
		if (_instance == NULL)
		{
			_instance = new ParticleFactory();
		}
		return _instance;
	}
	MovingParticle* makeParticle(CParticleF& p, MovingParticleType type, float tm)
	{
		MovingParticle* particle = new MovingParticle(p, type, tm);
		particles.push_back(particle);
		activeSet.insert(particle);
		return particle;
	}
	bool inactivate(MovingParticle* p)
	{
		set<MovingParticle*>::iterator it = activeSet.find(p);
		if (it == activeSet.end())
		{
			return false;
		}
		else
		{
			activeSet.erase(it);
			return true;
		}
	}
	void clean()
	{
		for (int i = 0; i<particles.size(); ++i)
		{
			delete particles[i];
		}
		particles.clear();
		MovingParticle::_id = 0;
	}
	MovingParticle* getNext()
	{
		float t = numeric_limits<float>::infinity();
		MovingParticle* p = NULL;
		for (set<MovingParticle*>::iterator it = activeSet.begin(); it != activeSet.end(); ++it)
		{
			MovingParticle* q = *it;
			if (q->event.t < t)
			{
				t = q->event.t;
				p = q;
			}
		}
		return p;
	}
	~ParticleFactory()
	{
		clean();
	}
	vector<MovingParticle*> particles;
	set<MovingParticle*> activeSet;
private:
	static ParticleFactory* _instance;
	ParticleFactory()
	{
		MovingParticle::_id = 0;
	}
};

/*
Utility functions
*/

bool
calculateBisectorVelocity(CParticleF a, CParticleF o, CParticleF b, float& vx, float& vy);

vector<MovingParticle*> vectorize(MovingParticle* p);

/*
true if the vector of points are in clock-wise direction
*/
bool clockWise(vector<MovingParticle*>& particles);

