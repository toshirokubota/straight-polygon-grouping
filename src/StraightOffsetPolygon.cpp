#include <vector>
#include <set>
using namespace std;

#include <StraightOffsetPolygon.h>
#include <IntersectionConvexPolygons.h>
#include <MovingParticle.h>
#include <Graph.h>
#include <szMiscOperations.h>

StraightOffsetPolygon::StraightOffsetPolygon(const vector<MovingParticle*>& vertices, float time)
{
	particles = vertices;
	this->created = time;
	//set the neighborhood
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
		MovingParticle::setNeighbors(p, r, q);
		pset.insert(p);
	}
	//calculate velocity if necessary
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		if (p->isInitialized()) continue;
		CParticleF o = p->getP();
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
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
		calculateBisectorVelocity(r->getP(), p->getP(), q->getP(), vx, vy);
		p->setVelocity(vx, vy);
	}
	//find the next event for each particle.
	//we only need to compute for those whose event has not been set of its participants no longer part of this polygon.
	for (int i = 0; i < particles.size(); ++i)
	{
		//MovingParticle* p = particles[i];
		//p->updateEvent();
		/*EventStruct event = p->getEvent();
		if (event.type == UnknownEvent ||
			pset.find(event.q) == pset.end() ||
			pset.find(event.r) == pset.end())
		{
			p->updateEvent(particles);
		}*/
	}
}

/*pair<StraightOffsetPolygon*, StraightOffsetPolygon*> 
StraightOffsetPolygon::applyNextEvent()
{
	//find the next event
	EventStruct ev(std::numeric_limits<float>::infinity(), UnknownEvent, NULL);
	for (int i = 0; i<particles.size(); ++i)
	{
		EventStruct ev0 = particles[i]->getEvent();
		if (ev0.t < ev.t)
		{
			ev = ev0;
		}
	}
	//move the particles
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->update(ev.t - particles[i]->getTime());
	}
	MovingParticle* p = (MovingParticle*)ev.p;
	pair<MovingParticle*, MovingParticle*> pp = p->applyEvent();
	pair<StraightOffsetPolygon*, StraightOffsetPolygon*> polys(NULL, NULL);
	if (pp.first)
	{
		polys.first = new StraightOffsetPolygon(MovingParticle::vectorize(pp.first), ev.t);
	}
	if (pp.second && 
		polys.first->exist(pp.second)==false //this particle is not a part of the first polygon.
		)
	{
		polys.second = new StraightOffsetPolygon(MovingParticle::vectorize(pp.second), ev.t);
	}
	return polys;
}*/

vector<CParticleF> 
StraightOffsetPolygon::snapShot() const
{
	vector<CParticleF> vp;
	for (int i = 0; i < particles.size(); ++i)
	{
		vp.push_back(particles[i]->getP());
	}
	return vp;
}

/*
This function checks for a self-intersecting side of a polygon
*/
bool
simplePolygonCheck(vector<MovingParticle*>& particles)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		CParticleF p = particles[i]->getP();
		CParticleF p2 = particles[i]->getNext()->getP();

		//NOTE: j starts j+2 to skip the adjacent side
		for (int j = i + 2; j < particles.size(); ++j)
		{
			CParticleF q = particles[j]->getP();
			CParticleF q2 = particles[j]->getNext()->getP();
			pair<float, float> param = _IntersectConvexPolygon::intersect(p, p2, q, q2);
			if (param.first>0 && param.first<1.0 && param.second>0 && param.second < 1.0)
			{
				return false;
			}
		}
	}
	return true;
}

/*
Trace a tree around it to construct particles representing vertices of a straight polygon.
The implementation is not pretty and may be simplified. But I am not sure at this time.
*/
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

/*
Trace a forrest and create a polygon for each tree.
*/
vector<StraightOffsetPolygon*>
traceForrest(vector<Vertex<CParticleF>*>& forrest)
{
	for (int i = 0; i < forrest.size(); ++i)
	{
		forrest[i]->Reset();
	}
	vector<StraightOffsetPolygon*> polygons;
	while (true)
	{

		vector<MovingParticle*> points;
		for (int j = 0; j < forrest.size(); ++j)
		{
			if (forrest[j]->color == White)
			{
				_trace(forrest[j], NULL, points);
				//remove a copy of the first one at the end of the trace
				//points.erase(points.end() - 1);
				StraightOffsetPolygon* poly = new StraightOffsetPolygon(points, 0.0f);
				polygons.push_back(poly);
				break;
			}
		}
		if (points.empty())
		{
			break;
		}
	}
	return polygons;
}

