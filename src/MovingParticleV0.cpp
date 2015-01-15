#include <MovingParticle.h>
#include <vector>
#include <set>
#include <limits>
using namespace std;
#include <mex.h> //need for premature exit when something goes wrong.

#include <szmexutilitytemplate.h>
#include <szMiscOperations.h>
#include <IntersectionConvexPolygons.h>

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
	}
	return tr;
}

bool
calculateBisectorVelocity(CParticleF a, CParticleF o, CParticleF b, float& vx, float& vy)
{
	if (o == a) //leaf node
	{
		float ang = GetVisualDirection(o.m_X, o.m_Y, b.m_X, b.m_Y) + PI/2.0;
		a.m_X = o.m_X + cos(ang);
		a.m_Y = o.m_Y + sin(ang);
	}
	else if (o == b) //leaf node
	{
		float ang = GetVisualDirection(o.m_X, o.m_Y, a.m_X, a.m_Y) - PI / 2.0;
		b.m_X = o.m_X + cos(ang);
		b.m_Y = o.m_Y + sin(ang);
	}
	CParticleF bs = bisector(o, a, b);
	float ang = GetVisualAngle2(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y);
	float cs = cos((PI - Abs(ang)) / 2.0);
	bool bret = true;
	if (cs < 0.01)
	{
		cs = 0.01;
		bret = false;
	}
	float len = 1.0f / cs;
	bs.m_X *= len;
	bs.m_Y *= len;
	vx = bs.m_X;
	vy = bs.m_Y;
	return bret;
}

bool clockWise(vector<MovingParticle*>& particles)
{
	vector<CParticleF> pnts;
	for (int i = 0; i<particles.size(); ++i)
	{
		pnts.push_back(particles[i]->getP());
	}
	return ClockWise(pnts);
}


/*
Find the time to hit for the 1st particle (P) to a slanted plane formed with Q and R.
Returns Inf if no intersection within the bounded polygon.
*/
float
MovingParticle::intersectPolygonAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r)
{
	CParticleF po = p->p;
	CParticleF qo = q->p;
	CParticleF ro = r->p;
	//perp to plane
	float eps = 1.0e-5;
	float t = std::numeric_limits<float>::infinity();
	{
		float d = 100.0f; //Max(Distance(po, qo), Distance(po, ro)); //find the bounding limit
		CParticleF poly[3];
		//make a 3D plane with a unit-slant
		poly[0] = qo; 
		poly[1] = ro;
		poly[2].m_X = qo.m_X + d * q->v[0];
		poly[2].m_Y = qo.m_Y + d * q->v[1];
		poly[2].m_Z = d;
		CParticleF z = perp2Plane(poly[0], poly[2], poly[3]);
		CParticleF u(p->v[0], p->v[1], 1.0f);
		float t0 = intersectPlaneAndLine(qo, z, po, u);
		if (t0 > 0)
		{
			//now check if the intersection is within the line segment
			pair<float, float> param = _IntersectConvexPolygon::intersect(p->p, p->move(t0), q->move(t0), r->move(t0));
			if (param.second > 0 && param.second <= 1.0f)
			{
				t = t0;
			}
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
	pair<float,float> param = _IntersectConvexPolygon::intersect(p0, p2, q0, q2);
	EventStruct ev(std::numeric_limits<float>::infinity(), CollisionEvent, this, this->next);
	if (param.first > 0)
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
		if (this->id == 1111 && q->id == 1109)
			ev.t += 0;
		if (this == q || this->prev == q) continue;
		//CParticleF w((q->v[0]+r->v[0])/2.0, (q->v[1]+r->v[1])/2.0); //vector perpendicular to the side
		//for the particle to hit the side, the particle cannot be moving away.
		//This condition tries to eliminate the two adjacent particles situated at the same place at the beginning of evolution.
		//float eta = -w.m_X * v[0] - w.m_Y * v[1];

		//if (eta <= 0) continue;

		float t0 = intersectPolygonAndLine(this, q, r); // period - p->time);
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
	event = findNextEdgeEvent();
	EventStruct ev2 = findNextSplitEvent();
	bool bret = false;
	/*if (ev1.t < event.t)
	{
		event = ev1;
		bret = true;
	}*/
	if (ev2.t < event.t)
	{
		event = ev2;
		bret = true;
	}
	return bret;
}

pair<MovingParticle*, MovingParticle*> 
MovingParticle::applyEvent()
{
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
		//setRelation(this, pnew);
		//setRelation(this->next, pnew);
	}
	else if (event.type == SplitEvent)
	{
		pnew[0] = factory->makeParticle(this->p, Split, event.t);
		pnew[1] = factory->makeParticle(this->p, Split, event.t);
		setNeighbors(pnew[0], event.p->prev, r);
		setNeighbors(pnew[1], q, event.p->next);
		factory->inactivate(p);
		///TK - need to implement this.
		///setRelation(p, pnew[0]);
		///setRelation(p, pnew[1]);
	}
	for (int i = 0; i<2; ++i)
	{
		if (pnew[i] == NULL) continue;
		float vx, vy;
		if (calculateBisectorVelocity(pnew[i]->prev->p, pnew[i]->p, pnew[i]->next->p, vx, vy) == false)
		{
			printf("failed to set velocity.\n");
			this->print();
			mexErrMsgTxt("applyEvent: failed to set velocity.");
		}
		pnew[i]->setVelocity(vx, vy);
	}
	return pair<MovingParticle*, MovingParticle*>(pnew[0], pnew[1]);
}
