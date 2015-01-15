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

void 
MovingParticle::print(char* tab) const
{
	printf("%s%d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f\n", tab, id, p0.m_X, p0.m_Y, p.m_X, p.m_Y, v[0], v[1]);
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

bool
calculateBisectorVelocity(CParticleF a, CParticleF o, CParticleF b, float& vx, float& vy)
{
	double eps = 1.0e-4;
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
MovingParticle::intersectSideAndLine(const MovingParticle* p, const MovingParticle* q, const MovingParticle* r)
{
	CParticleF po = p->p;
	CParticleF qo = q->p;
	CParticleF ro = r->p;
	//perp to plane
	float eps = 1.0e-4;
	float t = std::numeric_limits<float>::infinity();
	{
		if (isApproaching(po, p->v[0], p->v[1], qo, ro))
		{
			CParticleF qrv = perpDirection(qo, ro);
			float d = 100.0f; //Max(Distance(po, qo), Distance(po, ro)); //find the bounding limit
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

			CParticleF m;
			for (int i = 0; i < 4; ++i)
			{
				m.m_X += poly[i].m_X;
				m.m_Y += poly[i].m_Y;
				m.m_Z += poly[i].m_Z;
			}
			m.m_X /= 4.0f;
			m.m_Y /= 4.0f;
			m.m_Z /= 4.0f;

			CParticleF z = perp2Plane(m, poly[2], poly[3]);
			CParticleF u(p->v[0], p->v[1], 1.0f);
			float t0 = intersectPlaneAndLine(m, z, po, u);
			if (t0 > 0)
			{
				//now check if the intersection is within the line segment
				CParticleF p2 = p->move(t0);
				CParticleF q2a = q->move(t0);
				CParticleF r2a = r->move(t0);
				float dval = Distance2LineSegment(q2a, r2a, p2);
				if (dval < eps)
				{
					t = t0;
				}
				else
				{
					CParticleF q2b(qo.m_X + qrv.m_X*t0, qo.m_Y + qrv.m_Y*t0);
					CParticleF r2b(ro.m_X + qrv.m_X*t0, ro.m_Y + qrv.m_Y*t0);
					float dval2 = Distance2LineSegment(q2b, r2b, p2);
					if (dval2 < eps)
					{
						t = t0;
					}
				}
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
		if (this->id == 10 && (q->id == 43))
		{
			ev.t += 0;
		}
		if (this == q || this->prev == q) continue;

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
	if (event.type == CollisionEvent && next == event.q) return false;
	if (event.type == SplitEvent && event.q->bActive && event.q->next == event.r) return false;

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

bool 
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
		pnew[0]->anscestors[0] = this;
		pnew[0]->anscestors[1] = this->next;
	}
	else if (event.type == SplitEvent)
	{
		pnew[0] = factory->makeParticle(this->p, Split, event.t);
		pnew[1] = factory->makeParticle(this->p, Split, event.t);
		setNeighbors(pnew[0], event.p->prev, r);
		setNeighbors(pnew[1], q, event.p->next);
		factory->inactivate(p);
		///TK - need to implement this.
		pnew[0]->anscestors[0] = this;
		pnew[0]->anscestors[1] = NULL;
		pnew[1]->anscestors[0] = this;
		pnew[1]->anscestors[1] = NULL;
	}
	for (int i = 0; i<2; ++i)
	{
		if (pnew[i] == NULL) continue;
		float vx, vy;
		if (calculateBisectorVelocity(pnew[i]->prev->p, pnew[i]->p, pnew[i]->next->p, vx, vy) == false)
		{
			pnew[i]->setVelocity((q->v[0] + r->v[0]) / 2.0, (q->v[1] + r->v[1]) / 2);
			pnew[i]->bUnstable = true;
		}
		else
		{
			pnew[i]->setVelocity(vx, vy);
		}
	}
	return true;
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

vector<vector<CParticleF>> 
MovingParticle::takeSnapshots()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	set<MovingParticle*> pset;
	vector<vector<CParticleF>> snapshots;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = vectorize(p);
				vector<CParticleF> shot;
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
					shot.push_back(vp[i]->p);
				}
				bdone = false;
				snapshots.push_back(shot);
			}
		}
		if (bdone) break;
	}
	return snapshots;
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
					return false;
				}
			}
		}
	}
	return true;
}

void
_trace(Vertex<CParticleF>* u, Vertex<CParticleF>*prev, vector<MovingParticle*>& points)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* particle = factory->makeParticle(u->key, Initial, 0.0f);
	points.push_back(particle);
	u->color = Black;
	CParticleF pu = u->key;
	CParticleF qu = prev == NULL ? CParticleF(pu.m_X - 1.0, pu.m_Y) : prev->key;
	int count = 1;
	while (true)
	{
		Vertex<CParticleF>* v = NULL;
		float minAngle = 3 * PI;
		//Walk in clock-wise  order.
		for (int i = 0; i < u->aList.size(); ++i)
		{
			Vertex<CParticleF>* w = u->aList[i]->v;
			if (w->color == White)
			{
				float ang = GetVisualAngle2(qu.m_X, qu.m_Y, w->key.m_X, w->key.m_Y, pu.m_X, pu.m_Y);
				if (ang <= 0) ang += 2 * PI;
				if (ang < minAngle)
				{
					minAngle = ang;
					v = w;
				}
			}
		}
		if (v != NULL)
		{
			count++;
			_trace(v, u, points);
			MovingParticle* particle2 = factory->makeParticle(u->key, Initial, 0.0f);
			points.push_back(particle2); //every branch will result in a new particle.
		}
		else
		{
			break;
		}
		qu = v->key; //update where we are coming from
	}
	//at a leaf, we need to provide additional one to account for two corners
	if (count < 2)
	{
		MovingParticle* particle2 = factory->makeParticle(u->key, Initial, 0.0f);
		points.push_back(particle2);
	}
	//When this is the first node of the trace, the last one is extra unless this is also a leaf node.
	if (prev == NULL && count > 2)
	{
		MovingParticle* p = *(points.end() - 1);
		factory->inactivate(p);
		points.erase(points.end() - 1);
	}
}

bool
MovingParticle::initializePolygon(vector<MovingParticle*>& particles)
{
	//set the neighborhood
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
		MovingParticle::setNeighbors(p, r, q);
	}
	//calculate velocity if necessary
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		if (p->isInitialized()) continue;
		CParticleF o = p->p;
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
		CParticleF b = q->p;
		CParticleF a = r->p;
		if (o == a) //leaf node
		{
			float ang = GetVisualDirection(o.m_X, o.m_Y, b.m_X, b.m_Y) + PI / 2.0;
			a.m_X = o.m_X + cos(ang);
			a.m_Y = o.m_Y + sin(ang);
		}
		else if (o == b) //leaf node
		{
			float ang = GetVisualDirection(o.m_X, o.m_Y, a.m_X, a.m_Y) - PI / 2.0;
			b.m_X = o.m_X + cos(ang);
			b.m_Y = o.m_Y + sin(ang);
		}
		float vx, vy;
		calculateBisectorVelocity(a, o, b, vx, vy);
		p->setVelocity(vx, vy);
	}
	return true;
}

/*
Trace a forrest and create a polygon for each tree.
*/
bool
MovingParticle::traceForrest(vector<Vertex<CParticleF>*>& forrest)
{
	for (int i = 0; i < forrest.size(); ++i)
	{
		forrest[i]->Reset();
	}
	while (true)
	{

		vector<MovingParticle*> points;
		for (int j = 0; j < forrest.size(); ++j)
		{
			if (forrest[j]->color == White)
			{
				_trace(forrest[j], NULL, points);
				initializePolygon(points);
				break;
			}
		}
		if (points.empty())
		{
			break;
		}
	}
	return true;
}

/*void
MovingParticle::traceBack(MovingParticle* p, vector<MovingParticle*>& trace, set<MovingParticle*>& pset)
{
	if (pset.find(p) == pset.end())
	{
		pset.insert(p);
		if (p->type == Regular || p->type == Initial)
		{
			trace.push_back(p);
		}
		else
		{
			for (int i = 0; i < p->descendents.size(); ++i)
			{
				traceBack(p->descendents[i], trace, pset);
			}
		}
	}
}*/

/*
Given a polygon, it returns a series of initial particles by tracing particles backward in time.
*/
vector<MovingParticle*>
MovingParticle::traceBackPolygon(MovingParticle* p)
{
	vector<MovingParticle*> trace;
	set<MovingParticle*> pset;
	//traceBack(p, trace, pset);
	return trace;
}
