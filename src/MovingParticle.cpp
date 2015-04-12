#include <MovingParticle.h>
#include <vector>
#include <set>
#include <limits>
using namespace std;
#include <mex.h> //need for premature exit when something goes wrong.

#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <IntersectionConvexPolygons.h>
#include <Graph.h>

//enum MovingParticleType { Unknown, Initial, Regular, Merge, Collide, Split, Axis, Dummy };
MovingParticleType int2ParticleType(int i)
{
	MovingParticleType type = Unknown;
	switch (i)
	{
	case 1:
		type = Initial;
		break;
	case 2:
		type = Regular;
		break;
	case 3:
		type = Merge;
		break;
	case 4:
		type = Collide;
		break;
	case 5:
		type = Split;
		break;
	case 6:
		type = Axis;
		break;
	case 7:
		type = Dummy;
		break;
	}
	return type;
}

void 
MovingParticle::print(char* tab) const
{
	printf("%s%d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f\n", tab, id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1]);
}

vector<float>
MovingParticle::dump2vector()
{
	vector<float> v;
	v.push_back(id);
	v.push_back(p0.m_X);
	v.push_back(p0.m_Y);
	v.push_back(p.m_X);
	v.push_back(p.m_Y);
	v.push_back(this->v[0]);
	v.push_back(this->v[1]);
	v.push_back(prev == NULL ? -1 : prev->id);
	v.push_back(next == NULL ? -1 : next->id);
	v.push_back(type);
	v.push_back(created);
	v.push_back(time);
	v.push_back(reflexive);
	v.push_back(bActive ? 1.0f : 0.0f);
	v.push_back(bInitialized ? 1.0f : 0.0f);
	v.push_back(bUnstable ? 1.0f : 0.0f);
	v.push_back(parents[0] == NULL ? -1 : parents[0]->id);
	v.push_back(parents[1] == NULL ? -1 : parents[1]->id);
	v.push_back(children[0] == NULL ? -1 : children[0]->id);
	v.push_back(children[1] == NULL ? -1 : children[1]->id);
	v.push_back(event.type);
	v.push_back(event.q == NULL ? -1.0f : event.q->id);
	v.push_back(event.r == NULL ? -1.0f : event.r->id);
	v.push_back(event.t);
	return v;
}

/*
A reflex particle is non convex one that can split a side of a polygon.
*/
bool 
MovingParticle::isReflex() const //check if it is concave (reflexive) that allows splitting of a side;
{
	//return GetVisualAngle2(prev->p.m_X, prev->p.m_Y, next->p.m_X, next->p.m_Y, p.m_X, p.m_Y) <= 0;
	return this->reflexive <= 0.0f;
}

//find the angle of propagating front p-q.
float 
MovingParticle::frontPropAngle(MovingParticle* p, MovingParticle* q)
{
	CParticleF p0 = p->move(1);
	CParticleF p1 = p->p;
	CParticleF p2 = q->p;
	CParticleF p3 = q->move(1);

	return GetVisualAngle2(p0.m_X, p0.m_Y, p2.m_X, p2.m_Y, p1.m_X, p1.m_Y) + 
		GetVisualAngle2(p1.m_X, p1.m_Y, p3.m_X, p3.m_Y, p2.m_X, p2.m_Y);

}

vector<MovingParticle*>
MovingParticle::vectorize(MovingParticle* start)
{
	vector<MovingParticle*> tr;
	MovingParticle* p = start;
	set<MovingParticle*> pset;
	bool success = true;
	do
	{
		tr.push_back(p);
		if (pset.find(p) != pset.end()) //premature loop is found
		{
			success = false;
			break;
		}
		if (p->next == NULL)
		{
			success = false;
			break;
		}
		pset.insert(p);
		p = p->next;
	} while (p != start);
	if (success == false)
	{
		printf("failed vectorization data.\n");
		for (int i = 0; i<tr.size(); ++i)
		{
			printf("%d %f %f %d\n", i + 1, tr[i]->p.m_X, tr[i]->p.m_Y, tr[i]->id);
		}
		mexErrMsgTxt("vectorize: failed to vectorize a shape.");
		//tr.clear();
	}
	return tr;
}

//find the particle with the next event in line.
MovingParticle*
MovingParticle::getNextEvent()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* p = NULL;
	float time = std::numeric_limits<float>::infinity();
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); it++)
	{
		if ((*it)->event.t < time)
		{
			p = *it;
			time = p->event.t;
		}
	}
	return p;
}

/*
Hopefully a more stable way to compute the bisector velocity.
It uses a vector direction of a colliding line segment in (dx, dy) and 
a vector diction of a line segment being split in (ux, uy).
It then returns the velocity in (vx, vy).
*/
bool
calculateBisectorVelocity2(CParticleF o, float ox, float oy, float ux, float uy, float& vx, float& vy)
{
	CParticleF a(o.m_X + ox, o.m_Y + oy);
	CParticleF b(o.m_X + ux, o.m_Y + uy);
	CParticleF bs = bisector(o, a, b);
	double ang = GetVisualAngle2(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y);
	double cs = cos((PI - Abs(ang)) / 2.0);
	bool bret = true;
	if (cs < 0.01)
	{
		cs = 0.01;
		bret = false;
	}
	double len = 1.0f / cs;
	bs.m_X *= len;
	bs.m_Y *= len;
	vx = (float)bs.m_X;
	vy = (float)bs.m_Y;
	return bret;
}

bool
calculateBisectorVelocity(CParticleF a, CParticleF o, CParticleF b, float& vx, float& vy)
{
	double eps = 1.0e-3;
	if (Distance(o, a) < eps || Distance(o, b) < eps) //too close
	{
		return false;
	}
	CParticleF bs = bisector(o, a, b);
	double ang = GetVisualAngle2(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y);
	double cs = cos((PI - Abs(ang)) / 2.0);
	bool bret = true;
	if (cs < 0.01)
	{
		cs = 0.01;
		bret = false;
	}
	double len = 1.0f / cs;
	bs.m_X *= len;
	bs.m_Y *= len;
	vx = (float) bs.m_X;
	vy = (float) bs.m_Y;
	return bret;
}

int clockWise(vector<MovingParticle*>& particles)
{
	vector<CParticleF> pnts;
	for (int i = 0; i<particles.size(); ++i)
	{
		pnts.push_back(particles[i]->getP());
	}
	return ClockWise(pnts);
}

/*
This utility method check when a moving particle (P) with a moving line of Q and its next.
It returns the time of the split/collision.
*/
float 
MovingParticle::_splitTime(const MovingParticle* q, float eps) const
{
	MovingParticle* r = q->next;
	CParticleF po = this->p;
	CParticleF qo = q->p;
	CParticleF ro = r->p;

	CParticleF qrv = perpDirection(qo, ro);
	float d = 1.0f; //Max(Distance(po, qo), Distance(po, ro)); //find the bounding limit
	CParticleF poly[4];
	//make a 3D plane with a unit-slant
	poly[0] = qo;
	poly[1] = ro;
	poly[2].m_X = ro.m_X + d * qrv.m_X;
	poly[2].m_Y = ro.m_Y + d * qrv.m_Y;
	poly[2].m_Z = d;
	poly[3].m_X = qo.m_X + d * qrv.m_X;
	poly[3].m_Y = qo.m_Y + d * qrv.m_Y;
	poly[3].m_Z = d;

	float mind = std::numeric_limits<float>::infinity();
	int idx = -1;
	for (int i = 0; i < 4; ++i)
	{
		float dval = Distance(poly[i], po);
		if (dval < mind)
		{
			mind = dval;
			idx = i;
		}
	}

	CParticleF z1 = perp2Plane(poly[idx], poly[(idx + 1) % 4], poly[(idx + 2) % 4]);
	CParticleF z2 = perp2Plane(poly[idx], poly[(idx + 2) % 4], poly[(idx + 3) % 4]);
	CParticleF z3 = perp2Plane(poly[idx], poly[(idx + 3) % 4], poly[(idx + 1) % 4]);
	CParticleF z((z1.m_X + z2.m_X + z3.m_X) / 3, (z1.m_Y + z2.m_Y + z3.m_Y) / 3, (z1.m_Z + z2.m_Z + z3.m_Z) / 3);
	CParticleF u(this->v[0], this->v[1], 1.0f);
	float t0 = intersectPlaneAndLine(poly[idx], z, po, u);
	return t0;
}

/*
Check if the particle P will be located on the side formed by Q and its next at time t from the current time.
*/
bool
MovingParticle::_onSideAt(const MovingParticle* q, float t, float eps) const
{
	MovingParticle* r = q->next;
	CParticleF p2 = this->move(t);
	CParticleF q2a = q->move(t);
	CParticleF r2a = r->move(t);
	float dval = Distance2LineSegment(q2a, r2a, p2);
	return dval < eps;
}

/*
For a particle created by an event (CAUSE), set its parents. 
*/
void
MovingParticle::_setParents(EventStruct cause)
{
	MovingParticle* pe = (MovingParticle*)cause.p;
	MovingParticle* qe = (MovingParticle*)cause.q;
	MovingParticle* re = (MovingParticle*)cause.r;
	if (cause.type == CollisionEvent)
	{
		this->parents[0] = pe;
		this->parents[1] = qe;
	}
	else if (cause.type == SplitEvent)
	{
		this->parents[0] = pe;
		/*ParticleFactory* factory = ParticleFactory::getInstance();
		this->parents[0] = pe;
		MovingParticle* pe2 = factory->makeParticle(pe->p, Dummy, pe->time);
		pe2->event = pe->event;
		factory->inactivate(pe2);*/
		
		/*float dq = Distance(pe->p, qe->p);
		float dr = Distance(pe->p, re->p);
		MovingParticle* bridge = dq < dr ? qe : re;
		if (ClockWise(p.m_X, p.m_Y, pe->p.m_X, pe->p.m_Y, bridge->p.m_X, bridge->p.m_Y) > 0)
		{
			this->parents[0] = pe;
			this->parents[1] = bridge;
		}
		else
		{
			this->parents[0] = bridge;
			this->parents[1] = pe;
		}*/
	}
}

/*
Find the time to hit for the 1st particle (P) to a slanted plane formed with Q and R.
Returns Inf if no intersection within the bounded polygon.
*/
float
MovingParticle::intersectSideAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r)
{
	float t = std::numeric_limits<float>::infinity();
	float t0 = p->_splitTime(q);
	if (t0 > 0)
	{
		if (p->_onSideAt(q, t0))
		{
			t = t0;
		}
	}
	return t;
}

EventStruct
MovingParticle::findNextEdgeEvent() const
{
	CParticleF p2(p.m_X + v[0], p.m_Y + v[1]);
	CParticleF q0 = next->p;
	CParticleF q2(q0.m_X + next->v[0], q0.m_Y + next->v[1]);
	pair<float,float> param = _IntersectConvexPolygon::intersect(p, p2, q0, q2);
	EventStruct ev(std::numeric_limits<float>::infinity(), CollisionEvent, this, this->next);
	if (param.first > 0 && param.second > 0)
	{
		ev.t = time + param.first;
	}
	return ev;
}

EventStruct
MovingParticle::findNextSplitEvent() const
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	EventStruct ev(std::numeric_limits<float>::infinity(), SplitEvent, this);
	for (set<MovingParticle*>::iterator j = factory->activeSet.begin(); j != factory->activeSet.end(); ++j)
	{
		const MovingParticle* q = *j; 
		const MovingParticle* r = q->next;
		if (id == 431 && (q->id==430 || q->id==428))
			ev.t += 0;
		if (this == q || this->prev == q) continue;
		if ((Abs(this->created - q->created) < 1.0e-8 && Distance(this->p, q->p) < 1.0e-3) ||
			(Abs(this->created - r->created) < 1.0e-8 && Distance(this->p, r->p) < 1.0e-3))
		{
			//Nov. 15th, 2014
			//two particles created at the same time at the same place.
			//they cannot be splitting each other. However, because of numerical precision, we may get a very small positive 
			//collision time.
			continue;
		}

		//if (eta <= 0) continue;
		float t0 = intersectSideAndLine(this, q, r); // period - p->time);
		if (t0 >= 0)
		{
			float t = t0 + time;
			if (t < ev.t)
			{
				ev.t = t;
				ev.q = q;
				ev.r = r;
			}
		}
	}
	return ev;
}


bool
MovingParticle::updateEvent()
{
	//check if the event needs to be updated. If not, then simply return.
	if (this->bUnstable) return false; //unstable particle.
	if (this->bActive == false) return  false;
	bool bChanged = false;
	if (id == 82)
		id += 0;
	if (event.type == UnknownEvent)
	{
		if (isReflex())
		{
			event = findNextSplitEvent();
		}
		else
		{
			event = findNextEdgeEvent();
		}
		bChanged = true;
	}
	else if (event.type == CollisionEvent)
	{
		EventStruct event2 = findNextEdgeEvent();
		if (event.q != event2.q)
		{
			event = event2;
			bChanged = true;
		}
	}
	else if (event.type == SplitEvent)
	{
		if (event.q == NULL || event.r == NULL)
		{
			event = findNextSplitEvent();
			bChanged = true;
		}
		//Nov 25th.
		//The event time for split can get smaller in a situation like this.
		//Segment a-b is moving along the collision path, but the length gets smaller until either a or b experience collision.
		//At the time, the same segment's angle widen and it becomes possible to collide with c.
		//The speed of the side does not change but the visual angle can enlarge, thus enabling a collision.
		//else if (event.q->bActive && event.q->next == event.r)
		//{
			//no need to re-calculate
		//}
		else 
		{
			EventStruct event2 = findNextSplitEvent();
			if (event.q->bActive && event.q->next == event.r) //the current event is alive
			{
				if (event.t > event2.t) //only replace if it is strictly earlier
				{
					event = event2; 
					bChanged = true;
				}
			}
			else //the current event segment is destroyed, but the new one is a split sub-segment
			{
				bool b1 = false;
				float t1 = std::numeric_limits<float>::quiet_NaN();
				if (event.q->bActive)
				{
					t1 = _splitTime(event.q);
					if (t1 >= 0)
					{
						b1 = _onSideAt(event.q, t1);
					}
				}
				bool b2 = false;
				float t2 = std::numeric_limits<float>::quiet_NaN();
				if (event.r->bActive)
				{
					t2 = _splitTime(event.r->prev);
					if (t2 >= 0)
					{
						b2 = _onSideAt(event.r->prev, t2);
					}
				}
				if (b1)
				{
					event.r = event.q->next;
					event.t = t1 + time;
				}
				else if (b2)
				{
					event.q = event.r->prev;
					event.t = t2 + time;
				}
				else
				{
					event = event2;
				}
				bChanged = true;
			}
		}
	}
	if (bChanged == true)
	{
		if (event.q != NULL)
		{
			((MovingParticle*)event.q)->dependent.insert(this);
		}
		if (event.type == SplitEvent && event.r != NULL)
		{
			((MovingParticle*)event.r)->dependent.insert(this);
		}
	}
	return event.type != UnknownEvent;
}

#include <ParticleDirection.h>

bool 
MovingParticle::applyEvent()
{
	float eps = 1.0e-3;
	MovingParticle* p = this;
	MovingParticle* q = (MovingParticle*)event.q;
	MovingParticle* r = (MovingParticle*)event.r;
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* pnew[2] = { NULL, NULL };
	if (event.type == CollisionEvent)
	{
		pnew[0] = factory->makeParticle(p->p, Collide, event.t);
		setNeighbors(pnew[0], p->prev, p->next->next);
		factory->inactivate(p);
		factory->inactivate(q);
	}
	else if (event.type == SplitEvent)
	{
		pnew[0] = factory->makeParticle(this->p, Split, event.t);
		pnew[1] = factory->makeParticle(this->p, Split, event.t);

		setNeighbors(pnew[0], event.p->prev, r);
		setNeighbors(pnew[1], q, event.p->next);
		factory->inactivate(p);
	}

	for (int i = 0; i < 2; ++i)
	{
		if (pnew[i] == NULL) continue;

		pnew[i]->_setParents(event); //set parents of the new particle
		float ox = p

		float vx, vy;
		if (calculateBisectorVelocity2(pnew[i]->prev->p, pnew[i]->p, pnew[i]->next->p, vx, vy) == false)
		{
			pnew[i]->setVelocity(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
			//pnew[i]->setVelocity((q->v[0] + r->v[0]) / 2.0, (q->v[1] + r->v[1]) / 2);
			pnew[i]->bUnstable = true;
		}
		else
		{
			pnew[i]->setVelocity(vx, vy);
		}
	}

	//update dependency
	for (int i = 0; i < 2; ++i)
	{
		if (pnew[i] != NULL)
		{
			factory->updateQueue.insert(pnew[i]);
			factory->updateQueue.insert(pnew[i]->prev); //possibly a new collision into the new particle.
		}
	}
	if (event.type == CollisionEvent)
	{
		//if a new side has a wider angle, add other concave vertices to the update queue.
		if (frontPropAngle(this->prev, this) < frontPropAngle(pnew[0]->prev, pnew[0]) ||
			frontPropAngle(q, q->next) < frontPropAngle(pnew[0], pnew[0]->next))
		{
			for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
			{
				/////TK if ((*it)->isConcave())
				if ((*it)->event.type != CollisionEvent)
				{
					factory->updateQueue.insert(*it);
				}
			}
		}
	}
	for (set<MovingParticle*>::iterator it = dependent.begin(); it != dependent.end(); ++it)
	{
		factory->updateQueue.insert(*it);
	}
	for (set<MovingParticle*>::iterator it = q->dependent.begin(); it != q->dependent.end(); ++it)
	{
		factory->updateQueue.insert(*it);
	}
	((MovingParticle*)event.q)->removeDependent(this);
	if (event.type == SplitEvent && r != NULL)
	{
		for (set<MovingParticle*>::iterator it = r->dependent.begin(); it != r->dependent.end(); ++it)
		{
			factory->updateQueue.insert(*it);
		}
		((MovingParticle*)event.r)->removeDependent(this);
	}
	return true;
}

void 
MovingParticle::traceAndHandleUnstable(MovingParticle* p, vector<MovingParticle*>& unstable)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	if (p->bActive == false)
	{
		return;
	}
	if (p->prev == p->next)
	{
		return; //there are only two particles in this chain. quickFinish will take care of this.
	}
	float dp = Distance(p->p, p->prev->p);
	float dn = Distance(p->p, p->next->p);
	MovingParticle* q = NULL;
	if (dp < dn)
	{
		q = factory->makeParticle(p->prev->p, Merge, p->time);
		q->setNeighbors(q, p->prev->prev, p->next);
		q->_setParents(p->prev->event);
		factory->inactivate(p->prev);
	}
	else
	{
		q = factory->makeParticle(p->next->p, Merge, p->time);
		q->setNeighbors(q, p->prev, p->next->next);
		q->_setParents(p->next->event);
		factory->inactivate(p->next);
	}

	factory->updateQueue.insert(q);
	factory->updateQueue.insert(q->prev); //possibly a new collision into the new particle.
	for (set<MovingParticle*>::iterator it = p->dependent.begin(); it != p->dependent.end(); ++it)
	{
		factory->updateQueue.insert(*it);
	}
	for (set<MovingParticle*>::iterator it = p->prev->dependent.begin(); it != p->prev->dependent.end(); ++it)
	{
		factory->updateQueue.insert(*it);
	}
	for (set<MovingParticle*>::iterator it = p->next->dependent.begin(); it != p->next->dependent.end(); ++it)
	{
		factory->updateQueue.insert(*it);
	}

	if (q->initializeVelocity() == false)
	{
		traceAndHandleUnstable(q, unstable);
		factory->inactivate(q);
		//unstable.push_back(q);
	}
}

void
MovingParticle::removeUnstable()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	int count = 0;
	vector<MovingParticle*> unstable;
	while (true)
	{
		vector<MovingParticle*> Q;
		for (set<MovingParticle*>::iterator j = factory->activeSet.begin(); j != factory->activeSet.end(); ++j)
		{
			if ((*j)->bUnstable)
			{
				Q.push_back(*j);
			}
		}
		if (Q.empty())
		{
			break;
		}
		for (int i = 0; i < Q.size(); ++i)
		{
			count++;
			MovingParticle* p = Q[i];
			p->print("RM: ");
			traceAndHandleUnstable(p, unstable);
			factory->inactivate(p);
		}
	}
	if (count > 0)
	{
		//printf("RemoveUpstable: %d unstable particles handled.\n", count);
	}
}

bool
MovingParticle::initializeVelocity()
{
	bool bRet = false;
	if (bInitialized == false)
	{
		float vx, vy;
		if (calculateBisectorVelocity(this->prev->p, this->p, this->next->p, vx, vy) == false)
		{
			this->setVelocity(0.0f, 0.0f); // (thi->v[0] + r->v[0]) / 2.0, (q->v[1] + r->v[1]) / 2);
			this->bUnstable = true;
		}
		else
		{
			this->setVelocity(vx, vy);
			bRet = true;
		}
	}
	return bRet;
}

void
MovingParticle::quickFinish()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = vectorize(p);
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
				}
				if (vp.size() <= 3)
				{
					for (int i = 0; i < vp.size(); ++i)
					{
						factory->inactivate(vp[i]);
					}
				}
				bdone = false;
				break; //active set has changed. the interator needs to be initialized again.
			}
		}
		if (bdone) break;
	}
}


bool
MovingParticle::sanityCheck()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	bool bdone = true;
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		CParticleF p = (*it)->p;
		CParticleF p2 = (*it)->next->p;

		//NOTE: j starts j+2 to skip the adjacent side
		set<MovingParticle*>::iterator it2 = it;
		it2++;
		for (; it2 != factory->activeSet.end(); ++it2)
		{
			CParticleF q = (*it2)->p;
			CParticleF q2 = (*it2)->next->p;
			pair<float, float> param = _IntersectConvexPolygon::intersect(p, p2, q, q2);
			if (param.first>0.01 && param.first<0.99 && param.second>0.01 && param.second < 0.99)
			{
				//check if they are not coincidentally parallel to each other
				float dval = Distance2Line(q, q2, p);
				if(dval > 0.01)
				{

					printf("sanityCheck(): crossing at %d(%3.3f,%3.3f)-%d(%3.3f,%3.3f) and %d(%3.3f,%3.3f)-%d(%3.3f,%3.3f) with degree of [%3.3f,%3.3f] and %3.3f.\n",
						(*it)->id, p.m_X, p.m_Y, (*it)->next->id, p2.m_X, p2.m_Y,
						(*it2)->id, q.m_X, q.m_Y, (*it2)->next->id, q2.m_X, q2.m_Y,
						param.first, param.second, dval);
					(*it)->event.print();
					(*it)->next->event.print();
					(*it2)->event.print();
					(*it2)->next->event.print();
					return false;
				}
			}
		}
	}
	return true;
}

vector<vector<MovingParticle*>> 
MovingParticle::clusterParticles()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<vector<MovingParticle*>> clusters;
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = MovingParticle::vectorize(p);
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
				}
				bdone = false;
				clusters.push_back(vp);
			}
		}
		if (bdone) break;
	}
	return clusters;
}

void
_traceBack(MovingParticle* p, vector<MovingParticle*>& trace, set<MovingParticle*>& pset)
{
	pset.insert(p);
	if (p->getType() == Initial)
	{
		trace.push_back(p);
	}
	else
	{
		for (int i = 0; i < 2; ++i)
		{
			if (p->getParent(i) != NULL)
			{
				_traceBack(p->getParent(i), trace, pset);
			}
		}
	}
}

/*
Given a vectorized particles, it returns a series of initial particles by tracing particles backward in time.
*/
vector<MovingParticle*>
MovingParticle::traceBackPolygon(vector<MovingParticle*>& particles)
{
	vector<MovingParticle*> trace;
	set<MovingParticle*> pset;
	for (int i = 0; i<particles.size(); ++i)
	{
		if (pset.find(particles[i]) == pset.end())
		{
			_traceBack(particles[i], trace, pset);
		}
	}
	return trace;
}

void
_splitNow(vector<MovingParticle*>& points, //vectorized particles
		int beg, int end, //indices specifying the subsequence of points
		vector<int>& match, //if match[i]==j, then points[i] and points[j] overlaps
		vector<vector<MovingParticle*>>& polygons //add a region as being found.
		)
{
	if (beg >= end) return;
	vector<MovingParticle*> poly;
	for (int i = beg; i <= end;)
	{
		poly.push_back(points[i]);
		if (match[i] > i)
		{
			_splitNow(points, i + 1, match[i], match, polygons);
			i = match[i] + 1;
		}
		else
		{
			i++;
		}
	}
	if (poly.size() > 2)
	{
		polygons.push_back(poly);
	}
}

/*
From a vectorized particles, extract subsets that form closed regions.
*/
vector<vector<MovingParticle*>>
MovingParticle::closedRegions(vector<MovingParticle*>& points)
{
	vector<vector<MovingParticle*>> polygons;
	vector<int> match(points.size(), -1);
	float eps = 1.0e-5;
	for (int i = 0; i < points.size(); ++i)
	{
		int d = points.size();
		for (int j = 0; j < points.size(); ++j)
		{
			if (i == j) continue;
			if (Distance(points[i]->p0, points[j]->p0) < eps)
			{
				int k = Abs(i - j);
				if (k < d)
				{
					match[i] = j;
				}
			}
		}
	}

	_splitNow(points, 0, points.size() - 1, match, polygons);

	return polygons;
}
