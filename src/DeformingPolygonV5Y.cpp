#include <DeformingPolygonV5.h>
#include <map>
#include <set>
#include <szMiscOperations.h>
#include <szMexUtility.h>
#include <IntersectionConvexPolygons.h>
#include <DisjointSet.h>
#include <mex.h> //just for printf
#include <stack>

namespace StraightAxis
{

	EventStruct::EventStruct(float t, EventType type, MovingParticle* p, MovingParticle* q)
	{
		this->t = t;
		this->type = type;
		this->p = p;
		this->q = q;
		if(q == NULL)
		{
			r = NULL;
		}
		else
		{
			r = q->next;
		}
	}

/*
A unit vector from o to x in 2D
*/
CParticleF 
NormalizedDirection(CParticleF& x, CParticleF& o)
{
	float dx = x.m_X - o.m_X;
	float dy = x.m_Y - o.m_Y;
	float len = sqrt(dx*dx + dy*dy);
	if(len > 1.0e-5)
	{
		return CParticleF(dx/len, dy/len);
	}
	else
	{
		return CParticleF(0.0f, 0.0f);
	}
}

/*
A unit vector from o to x in 3D
*/
CParticleF 
NormalizedDirection3D(CParticleF& x, CParticleF& o)
{
	float dx = x.m_X - o.m_X;
	float dy = x.m_Y - o.m_Y;
	float dz = x.m_Z - o.m_Z;
	float len = sqrt(dx*dx + dy*dy + dz*dz);
	if(len > 1.0e-5)
	{
		return CParticleF(dx/len, dy/len, dz/len);
	}
	else
	{
		return CParticleF(0.0f, 0.0f, 0.0f);
	}
}

/*
return a unit vector that bisects the angle formed by three points: a-o-b.
*/
CParticleF
bisector(CParticleF& o, CParticleF& a, CParticleF& b)
{
	CParticleF x = NormalizedDirection(a, o);
	CParticleF y = NormalizedDirection(b, o);
	CParticleF z((x.m_X+y.m_X)/2, (x.m_Y+y.m_Y)/2);
	float vx = z.m_X;
	float vy = z.m_Y;
	float len0 = sqrt(vx*vx + vy*vy);
	if(len0 <= 1.0e-5) //this is a colinear point. 
	{
		float ang = GetVisualDirection(b.m_X, b.m_Y, a.m_X, a.m_Y) - PI/2.0;
		vx = cos(ang);
		vy = sin(ang);
	}
	else
	{
		vx = vx / len0;
		vy = vy / len0;
	}
	CParticleF bs(o.m_X+vx, o.m_Y+vy);
	//test for orientation
	vector<CParticleF> pnts(4);
	pnts[0] = o;
	pnts[1] = a;
	pnts[2] = bs;
	pnts[3] = b;
	if(ClockWise(pnts)<0)
	{
		vx = -vx;
		vy = -vy;
	}
	return CParticleF(vx, vy);
}

float
intersectPlaneAndLine(CParticleF op, CParticleF& z, CParticleF& ol, CParticleF u)
{
	//move the origin to ol.
	op.m_X -= ol.m_X;
	op.m_Y -= ol.m_Y;
	op.m_Z -= ol.m_Z;
	float t = (z.m_X * op.m_X + z.m_Y * op.m_Y + z.m_Z * op.m_Z) / (z.m_X * u.m_X + z.m_Y * u.m_Y + z.m_Z * u.m_Z);
	return t;
}

CParticleF
crossProduct(CParticleF& a, CParticleF& b)
{
	return CParticleF(a.m_Y*b.m_Z-a.m_Z*b.m_Y, a.m_Z*b.m_X-a.m_X*b.m_Z, a.m_X*b.m_Y-a.m_Y*b.m_X);
}

/*
Given three points in 3D, define a unit vector perpendicular to the plane going through the three points.
The vector is always upward (i.e. its Z-coordinate is always positive). But this should not really matter...
*/
CParticleF
perp2Plane(CParticleF& a, CParticleF& b, CParticleF& c)
{
	CParticleF u = NormalizedDirection3D(b, a);
	CParticleF v = NormalizedDirection3D(c, a);
	CParticleF z = crossProduct(u, v);
	if(z.m_Z < 0)
	{
		z.m_X = -z.m_X;
		z.m_Y = -z.m_Y;
		z.m_Z = -z.m_Z;
	}
	return z;
}

bool
calculateVelocity(CParticleF& a, CParticleF& o, CParticleF& b, float& vx, float& vy)
{
	CParticleF bs = bisector(o, a, b);
	float ang = GetVisualAngle2(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y);
	float cs = cos((PI - Abs(ang)) / 2.0);
	bool bret = true;
	if(cs < 0.1)
	{
		cs = 0.1;
		bret = false;
	}
	float len = 1.0f / cs;
	bs.m_X *= len;
	bs.m_Y *= len;
	vx = bs.m_X;
	vy = bs.m_Y;
	return bret;
}

bool 
MovingParticle::setVelocity()
{
	if(next && prev)
	{
		this->bStable = calculateVelocity(prev->p, p, next->p, v[0], v[1]);
	}
	else
	{
		mexErrMsgTxt("setVelocity: an isolated particle found.");		
		return false;
	}
}

void
setNeighbors(MovingParticle* p, MovingParticle* prev, MovingParticle* next)
{
	p->prev = prev;
	p->next = next;
	p->prev->next = p;
	p->next->prev = p;
}

void
setRelation(MovingParticle* child, MovingParticle* parent)
{
	if(child->type == Dummy)
	{
		if(Distance(child->p, child->prev->p) < Distance(child->p, child->next->p))
		{
			child = child->prev;
		}
		else
		{
			child = child->next;
		}
	}
	parent->descendents.push_back(child);
	child->anscestors.push_back(parent);
}

bool clockWise(vector<MovingParticle*>& particles)
{
	vector<CParticleF> pnts;
	for(int i=0; i<particles.size(); ++i)
	{
		pnts.push_back(particles[i]->p);
	}
	return ClockWise(pnts);
}

DeformingPolygon::DeformingPolygon(vector<CParticleF>& pnts, float tm, DeformingPolygon* pa)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	for(int i=0; i<pnts.size(); ++i)
	{
		MovingParticle* p = factory->makeParticle(pnts[i], Regular, tm);
	}
	_setMap(); //set particle-index map
	for(int i=0; i<particles.size(); ++i)
	{
		//particles[i]->p0 = particles[i]->p;
		int j = (i+1) % particles.size();
		particles[i]->next = particles[j];
		int k = (i-1+particles.size()) % particles.size();
		particles[i]->prev = particles[k];
		particles[i]->setVelocity();
	}
	created = tm;
	time = tm;
	orientation = ClockWise(pnts);
	id = _id++;
	parent = pa;
	died = 0;
}


DeformingPolygon::DeformingPolygon(vector<MovingParticle*>& pnts, float tm, DeformingPolygon* pa)
{
	this->particles = pnts;
	_setMap(); //set particle-index map
	for(int i=0; i<particles.size(); ++i)
	{
		//particles[i]->p0 = particles[i]->p;
		int j = (i+1) % particles.size();
		particles[i]->next = particles[j];
		int k = (i-1+particles.size()) % particles.size();
		particles[i]->prev = particles[k];
		particles[i]->setVelocity();
	}
	created = tm;
	time = tm;
	orientation = clockWise(particles);
	id = _id++;
	parent = pa;
	died = 0;
}

/*
Does this return true when the point is on the polygon?
*/
bool
inside2(const CParticleF& p, const vector<CParticleF>& pnts)
{
	bool c = false;
	float x = p.m_X;
	float y = p.m_Y;
	for (int i = 0, j = pnts.size()-1; i < pnts.size(); j = i++) 
	{
		float xi = pnts[i].m_X;
		float yi = pnts[i].m_Y;
		float xj = pnts[j].m_X;
		float yj = pnts[j].m_Y;
		if ( ((yi>y) != (yj>y)) &&
			(x <= (float)(xj-xi) * (float)(y-yi) / (float)(yj-yi) + xi) )
			c = !c;
	}
	return c;
}

EventStruct
NextEvent(MovingParticle* p, vector<MovingParticle*>& particles)
{
	float ptime = p->time;
	EventStruct ev(std::numeric_limits<float>::infinity(), UnknownEvent, p);
	{
		CParticleF ol = p->p;
		CParticleF ul(p->v[0], p->v[1], 1.0);
		for(int j=0; j<particles.size(); ++j)
		{
			if(p==particles[j] || p->prev==particles[j]) continue;
			if(p->id==1452 && particles[j]->id==1453)
				j+=0;
			CParticleF op = particles[j]->p;
			CParticleF oq = particles[j]->next->p;
			CParticleF or(op.m_X + particles[j]->v[0], op.m_Y + particles[j]->v[1], 1.0f);
			CParticleF up = perp2Plane(op, oq, or);
			float t0 = intersectPlaneAndLine(op, up, ol, ul);
			if(t0 >= 0) 
			{
				float t = t0 + ptime;
				if(t < ev.t)
				{
					float z = ol.m_Z + t0*ul.m_Z;
					CParticleF x(ol.m_X+t0*ul.m_X, ol.m_Y+t0*ul.m_Y);//the location of intersection in 2D
					CParticleF a(op.m_X + t0*particles[j]->v[0], op.m_Y + t0*particles[j]->v[1]);
					CParticleF b(oq.m_X + t0*particles[j]->next->v[0], oq.m_Y + t0*particles[j]->next->v[1]);
					float s = ((a.m_X-x.m_X)*(a.m_X-b.m_X) + (a.m_Y-x.m_Y)*(a.m_Y-b.m_Y))/((a.m_X-b.m_X)*(a.m_X-b.m_X) + (a.m_Y-b.m_Y)*(a.m_Y-b.m_Y));
					if(s>=0 && s<1.0)
					{
						if(particles[j] == p->next)
						{
							ev = EventStruct(t, CollisionEvent, p, particles[j]);
						}
						else
						{
							ev = EventStruct(t, SplitEvent, p, particles[j]);
						}
					}
				}
			}
		}
	}
	//ev.t = Max(ptime, ev.t);
	return ev;
}

int
DeformingPolygon::sanityCheck()
{
	int error_code = 0;
	map<MovingParticle*,bool> visited;
	{
		for(int i=0; i<particles.size(); ++i)
		{
			visited[particles[i]] = false;
		}
		int count = 0;
		MovingParticle* start = particles[0];
		MovingParticle* p = start;
		do
		{
			visited[p] = true;
			p = p->next;
			count++;
		}
		while(count < particles.size() && p != start);
		if(count < particles.size()) //premature ending
		{
			error_code = 1;
			return error_code;
		}
		else if(p != start) //unfinished
		{
			error_code = 2;
			return error_code;
		}
	}
	{
		for(int i=0; i<particles.size(); ++i)
		{
			visited[particles[i]] = false;
		}
		int count = 0;
		MovingParticle* start = particles[0];
		MovingParticle* p = start;
		do
		{
			visited[p] = true;
			p = p->prev;
			count++;
		}
		while(count < particles.size() && p != start);
		if(count < particles.size()) //premature ending
		{
			error_code = 3;
			return error_code;
		}
		else if(p != start) //unfinished
		{
			error_code = 4;
			return error_code;
		}
	}
	return error_code;
}

void
DeformingPolygon::quickFinish()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	if(particles.size()==1)
	{
		if(particles[0]->type != Axis)
		{
			MovingParticle* p = factory->makeParticle(particles[0]->p, Axis, time);
			setRelation(particles[0], p);
			particles[0]->bActive = false;
			p->bActive = false;
		}
	}
	if(particles.size()==2)
	{
		MovingParticleType type = particles[0]->type==Merge ? Merge: Axis;
		float x1 = (particles[0]->p.m_X + particles[1]->p.m_X) / 2;
		float y1 = (particles[0]->p.m_Y + particles[1]->p.m_Y) / 2;
		CParticleF q(x1, y1);
		MovingParticle* pa = factory->makeParticle(q, type, time);
		for(int i=0; i<2; ++i)
		{
			particles[i]->p = q;
			setRelation(particles[i], pa);
			particles[i]->bActive = false;
		}
		pa->bActive = false;
	}
	else if(particles.size()==3)
	{
		float t = std::numeric_limits<float>::quiet_NaN();
		MovingParticle* p = NULL;
		MovingParticle* q = NULL;
		MovingParticle* r = NULL;

		for(int n=0; n<3; ++n)
		{
			p = particles[n];
			q = particles[(n+1) % 3];
			r = particles[(n+2) % 3];
			pair<float,float> param = _IntersectConvexPolygon::intersect(p->p0, p->move(1.0), q->p0, q->move(1.0));
			if(Abs(param.first) < Distance(p->p, r->p) )
			{
				t = param.first; 
				break;
			}
		}
		MovingParticle* pnew = NULL;
		if(t==t){
			float x1 = (1.0-t) * p->p0.m_X + t * p->move(1.0).m_X;
			float y1 = (1.0-t) * p->p0.m_Y + t * p->move(1.0).m_Y;
			CParticleF f(x1,y1);
			pnew = factory->makeParticle(f, Axis, time);
		}
		else
		{
			CParticleF f((particles[0]->p.m_X+particles[1]->p.m_X+particles[2]->p.m_X)/3.0, (particles[0]->p.m_Y+particles[1]->p.m_Y+particles[2]->p.m_Y)/3.0);
			pnew = factory->makeParticle(f, Axis, time);
		}
		for(int i=0; i<3; ++i)
		{
			//particles[i]->p=f; ///TK
			setRelation(particles[i], pnew);
			particles[i]->bActive = false;
		}
		pnew->bActive = false;
	}
	died = time;
}

vector<MovingParticle*>
vectorize(MovingParticle* start)
{
	MovingParticle* p = start;
	vector<MovingParticle*> tr;
	set<MovingParticle*> pset;
	bool success = true;
	do
	{
		tr.push_back(p);
		if(pset.find(p) != pset.end()) //premature loop is found
		{
			success = false;
			break;
		}
		if(p->next == NULL)
		{
			success = false;
			break;
		}
		pset.insert(p);
		p = p->next;
	}
	while(p != start);

	if(success == false)
	{
		printf("failed vectorization data.\n");
		for(int i=0; i<tr.size(); ++i)
		{
			printf("%d %f %f %d\n", i+1, tr[i]->p.m_X, tr[i]->p.m_Y, tr[i]->id);
		}
		mexErrMsgTxt("vectorize: failed to vectorize a shape.");
	}
	return tr;
}

MovingParticle*
DeformingPolygon::applyCollision(MovingParticle* p, float delta)
{
	assert(pmap.find(p) != pmap.end() && pmap.find(p->next) != pmap.end());
	for(int i=0; i<particles.size(); ++i)
	{
		particles[i]->update(delta);
	}
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* pnew = factory->makeParticle(p->p, Collide, time + delta);
	setNeighbors(pnew, p->prev, p->next->next);
	pnew->setVelocity();
	setRelation(p, pnew);
	setRelation(p->next, pnew);
	time += delta;
	p->bActive = false;
	p->next->bActive = false;
	return pnew;
}

pair<MovingParticle*,MovingParticle*>
DeformingPolygon::applySplit(MovingParticle* p, MovingParticle* q, float delta)
{
	assert(pmap.find(p) != pmap.end() && pmap.find(q) != pmap.end());
	for(int i=0; i<particles.size(); ++i)
	{
		particles[i]->update(delta);
	}
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* pnew[2];
	pnew[0] = factory->makeParticle(p->p, Split, time + delta);
	pnew[1] = factory->makeParticle(p->p, Split, time + delta);
	setNeighbors(pnew[0], p->prev, q->next);
	setNeighbors(pnew[1], q, p->next);
	for(int i=0; i<2; ++i)
	{
		pnew[i]->setVelocity();
		setRelation(p, pnew[i]);
		setRelation(q, pnew[i]);
	}
	time += delta;
	p->bActive = false;
	return pair<MovingParticle*,MovingParticle*>(pnew[0], pnew[1]);
}

bool
affectedBy(EventStruct& current, EventStruct& another)
{
	vector<MovingParticle*> affected;
	if(current.type == CollisionEvent && another.type == CollisionEvent)
	{
		if(current.p==another.q || current.q==another.p)
		{
			return true;
		}
	}
	else if(current.type == CollisionEvent && another.type == SplitEvent)
	{
		if(current.p == another.q || current.p == another.r || current.q == another.p || current.q == another.q || current.q == another.r)
		{
			return true;
		}
	}
	else if(current.type == SplitEvent && another.type == CollisionEvent)
	{
		if(current.p == another.q  || current.q == another.p || current.q == another.q || current.r == another.p || current.r == another.q)
		{
			return true;
		}
	}
	else if(current.type == SplitEvent && another.type == SplitEvent)
	{
		if(current.p == another.q || current.p == another.r || 
			current.q == another.p || current.q == another.q || current.q == another.r ||
			current.r == another.p || current.r == another.q || current.r == another.r)
		{
			return true;
		}
	}
	return false;
}

bool need2Update(map<MovingParticle*,int>& pmap, EventStruct& ev)
{
	if(ev.type == CollisionEvent)
	{
		if(pmap.find(ev.p)!=pmap.end() && pmap.find(ev.q)!=pmap.end())
		{
			return false;
		}
		else return true;
	}
	else if(ev.type == SplitEvent)
	{
		if(pmap.find(ev.p)!=pmap.end() && pmap.find(ev.q)!=pmap.end() && pmap.find(ev.r)!=pmap.end())
		{
			return false;
		}
		else return true;
	}
	else
	{
		return true;
	}
}

/*
remove duplicates and kinks. 
*/
void
DeformingPolygon::removeUnstables()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	while(particles.size()>=3)
	{
		bool bChanged = false;
		for(int i=0; i<particles.size(); ++i)
		{
			if(particles[i]->bStable == false)
			{
				MovingParticle* p0 = particles[i];
				MovingParticle* q = NULL;
				MovingParticle* p = NULL;
				if(Distance(p0->p, p0->prev->p) < Distance(p0->p, p0->next->p))
				{
					p = p0->prev;
					q = factory->makeParticle(p->p, p->type, time);
					setNeighbors(q, p->prev, p0->next);
				}
				else
				{
					p = p0->next;
					q = factory->makeParticle(p->p, p->type, time);
					setNeighbors(q, p0->prev, p->next);
				}
				p0->bActive = false;
				p->bActive = false;
				setRelation(p, q);
				setRelation(particles[i], q);
				q->setVelocity();
				bChanged = true;
				particles = vectorize(q);

				break;
			}
		}
		if(bChanged == false) break;
	}
	for(int i=0; i<particles.size(); ++i)
	{
		if(particles[i]->bActive == false)
		{
			particles.erase(particles.begin() + i);
		}
	}

	_setMap();
}

vector<DeformingPolygon*> 
DeformingPolygon::deform(float period)
{
	int poid = 797; //73; //836; //122; //135;
	int paid = 1452;
	if(id==poid)
		poid += 0;
	vector<DeformingPolygon*> polygons(1, this);
	if(particles.size()<=3)
	{
		takeSnapshot();
		quickFinish();
		return polygons;
	}
	else
	{
		removeUnstables();
		takeSnapshot();
	}
	vector<MovingParticle*> queue;
	for(int i=0; i<particles.size(); ++i)
	{
		if(need2Update(pmap, particles[i]->event))
		{
			particles[i]->event = NextEvent(particles[i], particles);
			//particles[i]->event.t += time; //make it the polygon time
		}
		queue.push_back(particles[i]);
		if(particles[i]->id == paid)
		{
			MovingParticle* p = particles[i];
			EventStruct ev = p->event;
			EventStruct ev2 = NextEvent(particles[i], particles);
			i+=0;
			if(id==poid)
				i+=0;
			EventStruct ev3 = NextEvent(particles[i], particles);
		}
	}
	sort(queue.begin(), queue.end(), compareByEvent);
	for(int i=0; i<queue.size(); ++i)
	{
		MovingParticle* particle = queue[i];
		if(particle->bActive == false) continue;
		EventStruct current = queue[i]->event;
		if(current.t > period) break;
		float delta = current.t - time;
		if(particle->id == paid)
			i+=0;
		if(current.type == CollisionEvent)
		{
			MovingParticle* newp = applyCollision(current.p, delta);
			//update particles
			particles = vectorize(newp);
			_setMap();
			if(newp->bStable==false)
			{
				removeUnstables();
			}
			else
			{
				newp->event = NextEvent(newp, particles);
				//newp->event.t += time;
				queue.push_back(newp);
			}
			for(int j=queue.size()-1; j>i; --j)
			{
				//if(affectedBy(current, queue[j]->event))
				if(need2Update(pmap, queue[j]->event))
				{
					if(queue[j]->id == paid)
						paid += 0;
					EventStruct ev(std::numeric_limits<float>::infinity());
					if(queue[j]->bActive) 
					{
						ev = NextEvent(queue[j], particles);
						//ev.t += time;
						queue[j]->event = ev;
					}
					else
					{
						queue.erase(queue.begin() + j);
					}
				}
			}
			sort(queue.begin()+i+1, queue.end(), compareByEvent);
			if(particles.size()<=3)
			{
				quickFinish();
				break;
			}
		}
		else //split
		{
			pair<MovingParticle*,MovingParticle*> newpair = applySplit(current.p, current.q, delta);
			MovingParticle* pnew[2] = {NULL, NULL};
			pnew[0] = newpair.first;
			pnew[1] = newpair.second;
			assert(pnew[0] && pnew[1]);
			for(int k=0; k<2; ++k)
			{
				vector<MovingParticle*> ps = vectorize(pnew[k]);
				DeformingPolygon* newpoly = new DeformingPolygon(time, clockWise(ps), this);
				newpoly->particles = ps;
				newpoly->_setMap();
				children.push_back(newpoly);
				vector<DeformingPolygon*> newpolygons = newpoly->deform(period);
				polygons.insert(polygons.end(), newpolygons.begin(), newpolygons.end());
			}
			break;
		}
	}
	for(int i=0; i<polygons.size(); ++i)
	{
		time = Max(time, polygons[i]->time);
	}
	died = time;
	return polygons;
}

void
_trace(Vertex<CParticleF>* u, Vertex<CParticleF>*prev, vector<MovingParticle*>& points)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* particle = factory->makeParticle(u->key, Initial, 0.0f);
	points.push_back(particle);
	u->color = Black;
	CParticleF pu = u->key;
	CParticleF qu = prev==NULL ? CParticleF(pu.m_X-1.0, pu.m_Y): prev->key;
	while(true)
	{
		Vertex<CParticleF>* v = NULL;
		float minAngle = 2*PI;
		for(int i=0; i<u->aList.size(); ++i)
		{
			Vertex<CParticleF>* w = u->aList[i]->v;
			if(w->color == White)
			{
				float ang = GetVisualAngle2(qu.m_X, qu.m_Y, w->key.m_X, w->key.m_Y, pu.m_X, pu.m_Y);
				if(ang <= 0) ang += 2*PI;
				if(ang < minAngle) 
				{
					minAngle = ang;
					v = w;
				}
			}
		}
		if(v != NULL)
		{
			_trace(v, u, points);
			points.push_back(particle);
		}
		else
		{
			break;
		}
		qu = v->key;
	}
}

vector<MovingParticle*>
traceTree(vector<Vertex<CParticleF>*>& tree)
{
	for(int i=0; i<tree.size(); ++i)
	{
		tree[i]->Reset();
	}

	vector<MovingParticle*> points;
	_trace(tree[0], NULL, points);
	//remove a copy of the first one at the end of the trace
	points.erase(points.end()-1);
	for(int i=0; i<points.size(); ++i)
	{
		int i0 = (i-1+points.size()) % points.size();
		int i2 = (i+1) % points.size();
		setNeighbors(points[i], points[i0], points[i2]);
	}
	return points;
}

vector<MovingParticle*>
offsetPolygon(vector<MovingParticle*>& points, float margin, float eps)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<MovingParticle*> poly;
	for(int i=0; i<points.size(); ++i)
	{
		int i0 = (i-1+points.size()) % points.size();
		int i2 = (i+1) % points.size();
		CParticleF o = points[i]->p;
		CParticleF a = points[i0]->p;
		CParticleF b = points[i2]->p;
		if(Distance(a, o) < eps || Distance(b, o) < eps) 
		{
			continue;
		}
		if(Distance(a, b) > eps)
		{
			float vx, vy;
			calculateVelocity(a, o, b, vx, vy);
			//CParticleF q(o.m_X + margin * vx, o.m_Y + margin * vy);
			MovingParticle* particle = factory->makeParticle(o, Regular, 0.0f); //points[i];
			//particle->p0 = o;
			particle->v[0] = vx; 
			particle->v[1] = vy;
			setRelation(points[i], particle);
			poly.push_back(particle);
		}
		else //this is a leaf node of the tree.
		{
			float ang = GetVisualDirection(o.m_X, o.m_Y, a.m_X,  a.m_Y);
			float angles[] = {ang+PI/4.0, ang-PI/4.0};
			MovingParticle* q[2];
			q[0] = points[i];
			q[1] = factory->makeParticle(o, Regular, 0.0f);
			for(int j=0; j<2; ++j)
			{
				//float vx = sqrt(2.0) * margin*cos(angles[j]);
				//float vy = sqrt(2.0) * margin*sin(angles[j]);
				float vx = sqrt(2.0) * cos(angles[j]);
				float vy = sqrt(2.0) * sin(angles[j]);
				//CParticleF q(o.m_X+vx, o.m_Y+vy);
				MovingParticle* particle = factory->makeParticle(o, Regular, 0.0f);
				//particle->p0 = o;
				particle->v[0] = vx;
				particle->v[1] = vy;
				setRelation(points[i], particle);
				poly.push_back(particle);
			}
		}
	}
	for(int i=0; i<poly.size(); ++i)
	{
		//setting the neighbors, except velocity
		int j = (i+1) % poly.size();
		int k = (i-1+poly.size()) % poly.size();
		MovingParticle* p = poly[i];
		MovingParticle* prev = poly[k];
		MovingParticle* next = poly[j];
		setNeighbors(p, prev, next);
	}
	for(int i=0; i<poly.size(); ++i)
	{
		poly[i]->update(margin);
	}
	return poly;
}

mxArray*
StoreParticles(const vector<vector<MovingParticle*>>& polygons)
{
	const int dims[] = {polygons.size()};
	mxArray* cell = mxCreateCellArray(1, (mwSize*) dims);
	for(int i=0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size() ;
		const int dimsC[] = {n, 6};
		mxArray* ar = mxCreateNumericArray(2, (mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = polygons[i][j]->p0.m_X;
			p[n+j] = polygons[i][j]->p0.m_Y;
			p[2*n+j] = polygons[i][j]->p.m_X;
			p[3*n+j] = polygons[i][j]->p.m_Y;
			p[4*n+j] = polygons[i][j]->id;
			p[5*n+j] = polygons[i][j]->type;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

}