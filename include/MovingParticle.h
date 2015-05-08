#pragma once
#include <vector>
#include <set>
#include <map>
using namespace std;

#include <szParticleF.h>
#include <EventStruct.h>
#include <MiscGeometry.h>
#include <szMiscOperations.h>
#include <Graph.h>
#include <ParticleDirection.h>

enum MovingParticleType { Unknown, Initial, Regular, Merge, Collide, Split, Axis, Dummy };
struct ParticleFactory;
class ParticleSimulator;

MovingParticleType int2ParticleType(int i);
const int ParticleDumpSize = 23;

class MovingParticle
{
public:
	friend ParticleFactory; //so that the factory can access _id.
	friend ParticleSimulator;

	/*bool setVelocity(float vx, float vy)
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
	}*/
	void update(float delta)
	{
		if (bActive && !bUnstable)
		{
			p = move(delta);
			time += delta;
		}
	}
	void print(char* tab = "") const;

	CParticleF getP0() const { return p0; }
	CParticleF getP() const { return p; }
	void getVelocity(float& vx, float& vy) { vx = v[0]; vy = v[1]; }
	MovingParticle* getNext() const { return next; }
	MovingParticle* getPrev() const { return prev; }
	int getId() const { return id; }
	EventStruct getEvent() const { return event; }
	float getTime() const { return time; }
	float getCreatedTime() const { return created; }
	MovingParticleType getType() const { return type; }
	MovingParticle* getChildren(int j) const
	{
		if (j == 0) return this->children[0];
		else if (j == 1) return this->children[1];
		else return NULL;
	}
	MovingParticle* getParent(int j) const
	{
		if (j == 0) return this->parents[0];
		else if (j == 1) return this->parents[1];
		else return NULL;
	}
	bool isActive() const { return bActive; }
	//bool isInitialized() const { return bInitialized; }
	bool isReflex() const; //check if it is concave (reflex) that allows splitting of a side
	CParticleF move(float delta) const
	{
		return CParticleF(p.m_X + delta*v[0], p.m_Y + delta*v[1]);
	}
	CParticleF project(float time) const
	{
		float delta = time - created;
		return CParticleF(p0.m_X + delta*v[0], p0.m_Y + delta*v[1]);
	}
	bool updateEvent();
	//bool initializeVelocity();
	bool applyEvent();
	vector<float> dump2vector(); //store the current state as a vector of float 

	static float frontPropAngle(MovingParticle* p, MovingParticle* q); //find the angle of propagating front p-q.
	static void setNeighbors(MovingParticle* p, MovingParticle* prev, MovingParticle* next)
	{
		p->prev = prev;
		p->next = next;
		p->prev->next = p;
		p->next->prev = p;
		p->reflexive = GetVisualAngle2(p->prev->p.m_X, p->prev->p.m_Y, p->next->p.m_X, p->next->p.m_Y, p->p.m_X, p->p.m_Y);
		p->dependent.insert(prev);
		p->dependent.insert(next);
		prev->dependent.insert(p);
		next->dependent.insert(p);
		p->rear = prev->front;
		p->front = next->rear;
	}
	static vector<MovingParticle*> vectorize(MovingParticle* p);
	static MovingParticle* getNextEvent();
	static vector<MovingParticle*> getNextEvents(float eps = 1.0e-6);
	//static vector<vector<CParticleF>> takeSnapshots();
	static void removeUnstable(); //remove kinks and duplicates.
	static void quickFinish(); //remove polygons with less than 4 vertices.
	static bool sanityCheck();
	static vector<vector<MovingParticle*>> clusterParticles();
	//static void _traceBack(MovingParticle* p, vector<MovingParticle*>& trace, set<MovingParticle*>& pset);
	static vector<MovingParticle*> MovingParticle::traceBackPolygon(vector<MovingParticle*>& particles);
	static vector<vector<MovingParticle*>> closedRegions(vector<MovingParticle*>& points);
	static pair<MovingParticle*, float> findIntersection(MovingParticle* p, MovingParticle* q);

private:
	MovingParticle(CParticleF& p = CParticleF(0, 0, 0), MovingParticleType t = Unknown, float tm = 0.0f)
	{
		this->p = p;
		p0 = p;
		v[0] = v[1] = std::numeric_limits<float>::quiet_NaN();
		next = NULL;
		prev = NULL;
		parents[0] = parents[1] = NULL;
		children[0] = children[1] = NULL;
		id = _id++;
		type = t;
		created = tm;
		time = created;
		bActive = true;
		//bInitialized = false;
		bUnstable = false;
		//bNeedUpdate = true;
		init_event_time = -1;
		reflexive = std::numeric_limits<float>::quiet_NaN();
		event = EventStruct(std::numeric_limits<float>::infinity(), UnknownEvent, this);
	}

	float MovingParticle::_splitTime(const MovingParticle* q, const MovingParticle* r, float eps = 1.0e-3) const;
	bool _onSideAt(const MovingParticle* q, float t, float eps = 1.0e-3) const;
	void _setParents(EventStruct cause);
	//bool calculateVelocity(bool leaf = false);
	bool calculateVelocityR(float vp=1.0f, float vn=1.0f);
	bool initializeVelocity();

	static float intersectSideAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r);
	static void traceAndHandleUnstable(MovingParticle* p, vector<MovingParticle*>& removed);
	static bool correctOvershoot(MovingParticle* p, MovingParticle* q, pair<float, float> param);
	EventStruct findNextEdgeEvent() const;
	EventStruct findNextSplitEvent() const;
	bool addDependent(MovingParticle* p)
	{
		dependent.insert(p);
	}
	bool removeDependent(MovingParticle* p)
	{
		set<MovingParticle*>::iterator it = dependent.find(p);
		if (it == dependent.end())	return false;
		else
		{
			dependent.erase(it);
			return true;
		}
	}

	CParticleF p0;
	CParticleF p;
	MovingParticle* next;
	MovingParticle* prev;
	MovingParticle* parents[2];
	MovingParticle* children[2];
	ParticleDirection rear;
	ParticleDirection front;

	int id;
	MovingParticleType type;
	float created; //time this is created.
	float time;
	float init_event_time;
	//bool bInitialized; //false until its velocity is set.
	bool bActive; //true if it is still moving.
	bool bUnstable; //true if the bisector angle was too tight for accurate velocity calculation.
	//bool bNeedUpdate; //true if the event may need to be udpated.
	float reflexive; //measure of reflexivenessbReflexive
	EventStruct event;
	float v[2];
	set<MovingParticle*> dependent;
	static int _id;
};

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
		pmap[particle->id] = particle;
		if (particle->id == 443)
		{
			particle->id += 0;
		}
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

			(*it)->bActive = false;
			if ((*it)->id == 487)
				p->id += 0;
			activeSet.erase(it);
			return true;
		}
	}
	//when initialized from snapshots, some future particles have to be physically removed from memory.
	bool remove(MovingParticle* p)
	{
		particles.erase(find(particles.begin(), particles.end(), p));
		pmap.erase(p->id);
		set<MovingParticle*>::iterator it = find(activeSet.begin(), activeSet.end(), p);
		if (it != activeSet.end())
		{
			activeSet.erase(it);
		}
		delete p;
		return true;
	}
	void clean()
	{
		activeSet.clear();
		for (int i = 0; i<particles.size(); ++i)
		{
			delete particles[i];
		}
		particles.clear();
		pmap.clear();
		MovingParticle::_id = 0;
	}
	MovingParticle* get(int id)
	{
		return pmap[id];
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
	map<int, MovingParticle*> pmap;
	set<MovingParticle*> updateQueue;
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
int clockWise(vector<MovingParticle*>& particles);

