#include <DeformingPolygonV3.h>
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

/*
return true if x is on a line segment between [p, q) or (x0, x) crosses [p, q) where x0 is the original position of x.
It returns the point on the segment or the intersection in y. y is meaningful only when the return value is true.
*/
bool
overlap(CParticleF& x, CParticleF& p, CParticleF& q, CParticleF& x0, CParticleF& y, float eps)
{
	float a = p.m_X;
	float b = q.m_X-p.m_X;
	float c = p.m_Y;
	float d = q.m_Y-p.m_Y;
	float det = b*b + d*d;
	if(det < 1.0e-5) //p and q are too close
	{
		if(Distance(x, p) < eps)
		{
			y = x;
			return true;
		}
		else
			return false;
	}
	/*float t = (-a*b - c*d + b*x.m_X + d*x.m_Y)/(b*b + d*d);
	CParticleF y0(a+b*t, c+d*t); //a point on the line p-q that is closest to x
	if(t >=0 && t<1.0)
	{
		if(Distance(x, y0) < eps)
		{
			y = y0;
			return true;
		}
	}*/

	//check if (x0, x) crossses (p,q)
	pair<float,float> p1 = _IntersectConvexPolygon::intersect(x0, x, p, q);
	if(p1.first > 0 && p1.second>=0 && p1.second<=1.0) 
	{
		//(1-s)*p1 + s*p2
		float s = p1.first;
		y = CParticleF((1.0-s)*x0.m_X+s*x.m_X, (1.0-s)*x0.m_Y+s*x.m_Y);
		if(p1.first < 1.0 || Distance(x, y) < eps)
		{
			return true;
		}
	}
	//check if (x2, x) crossses (p,q)
	//if all failed
	return false;
}

bool
coLinear(CParticleF& a, CParticleF& b, CParticleF& c, float precision = 1.0e-3)
{
	float d = Distance2Line(a, c, b);
	return d < precision;
}

bool
coLinear(CParticleF& a, CParticleF& b, CParticleF& c, CParticleF& d, float precision = 1.0e-3)
{
	float d1 = Distance2Line(a, c, b);
	float d2 = Distance2Line(a, c, d);
	return d1 < precision && d2 < precision;
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
	float len = 1.0f / cs;
	bs.m_X *= len;
	bs.m_Y *= len;
	/*if(ang<0)
	{
	bs.m_X = -bs.m_X;
	bs.m_Y = -bs.m_Y;
	}*/
	vx = bs.m_X;
	vy = bs.m_Y;
	if(Abs(vx)>100000 || Abs(vy)>100000)
		vx += 0;
	return true;
}

bool
collided(MovingParticle* p, 
	float eps = 0.001 //minimum cosine that can be used to reliably compute the velocity.
	)
{
	if(p->next && p->prev)
	{
		float ang = GetVisualAngle2(p->next->p.m_X, p->next->p.m_Y, p->prev->p.m_X, p->prev->p.m_Y, p->p.m_X, p->p.m_Y);
		float cs = cos((PI - Abs(ang)) / 2.0);
		if(Abs(cs) < eps) return true;
		else return false;
	}
	else
	{
		return true;
	}
}

bool 
MovingParticle::setVelocity()
{
	float eps = 1.0e-5;
	if(next && prev)
	{
		/*if(Distance(next->p, prev->p) < eps)
		{
		CParticleF d = NormalizedDirection(next->p, p);
		v[0] = d.m_X;
		v[1] = d.m_Y;
		}
		else
		{
		calculateVelocity(prev->p, p, next->p, v[0], v[1]);
		}*/
		if(collided(this)) 
		{
			v[0] = v[1] = 0.0f;
			return false;
		}
		else
		{
			calculateVelocity(prev->p, p, next->p, v[0], v[1]);
			return true;
		}
	}
	else
	{
		v[0] = 0; v[1] = 0;
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
	//p->setVelocity();
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

DeformingPolygon::DeformingPolygon(vector<CParticleF>& pnts, float tm, int pa_id)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	float eps = 1.0e-5;
	for(int i=0; i<pnts.size(); ++i)
	{
		MovingParticle* p = factory->makeParticle(pnts[i], Regular, tm);
		particles.push_back(p);
	}
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
	id = _id++;
	parent_id = pa_id;
}

pair<float,MovingParticle*>
DeformingPolygon::nextEdgeEvent()
{
	float t = std::numeric_limits<float>::infinity();
	pair<float,MovingParticle*> res(t, (MovingParticle*)NULL);
	for(int i=0; i<particles.size(); ++i)
	{
		CParticleF p2 = particles[i]->move(1);
		CParticleF q2 = particles[i]->next->move(1);
		pair<float,float> param = _IntersectConvexPolygon::intersect(particles[i]->p, p2, particles[i]->next->p, q2);
		float t0 = Max(param.first, param.second);
		if(param.first>0 && param.second >0 && t0 < t)
		{
			t = t0;
			res.first = t;
			res.second = particles[i];
		}
	}
	return res;
}

pair<float,MovingParticle*>
DeformingPolygon::nextSplitEvent()
{
	float t = std::numeric_limits<float>::infinity();
	pair<float,MovingParticle*> res(t, (MovingParticle*)NULL);
	for(int i=0; i<particles.size(); ++i)
	{
		//CParticleF bs = bisector(particles[i]->p, particles[i]->prev->p, particles[i]->next->p);
		CParticleF dn = NormalizedDirection(particles[i]->next->p, particles[i]->p);
		CParticleF dp = NormalizedDirection(particles[i]->prev->p, particles[i]->p);
		if(dn.m_X * particles[i]->v[0] + dn.m_Y * particles[i]->v[1] < 0 &&
			dp.m_X * particles[i]->v[0] + dp.m_Y * particles[i]->v[1] < 0) //concave
		{
			CParticleF ol = particles[i]->p;
			CParticleF ul(particles[i]->v[0], particles[i]->v[1], 1.0);
			for(int j=0; j<particles.size(); ++j)
			{
				if(particles[i]==particles[j] || particles[i]->next==particles[j] || particles[i]->prev==particles[j]) continue;
				CParticleF op = particles[j]->p;
				CParticleF oq = particles[j]->next->p;
				CParticleF or(op.m_X + particles[j]->v[0], op.m_Y + particles[j]->v[1], 1.0f);
				CParticleF up = perp2Plane(op, oq, or);
				float t0 = intersectPlaneAndLine(op, up, ol, ul);
				if(t0 > 0 && t0<t) 
				{
					/*check if the intersection point is within the sweep of the line segment.
					For the sweeping line segment p-q, the intersection point (x) has to be on the same side of the moving line with 
					the other end point. It leads to four conditions: 
					1) the visual angle of x-p-q has to be smaller than the visual angle of p'-p-q where p' is a point on the p's projectile, 
					2) the visual angle of x-q-p has to be smaller than the visual angle of q'-q-p where q' is a point on the q's projectile,
					3) x is on the same side with p' with respect to p-q, and
					4) x is on the same side with q' with respect to p-q.*/
					CParticleF x(ol.m_X+t0*ul.m_X, ol.m_Y+t0*ul.m_Y);//the location of intersection in 2D
					CParticleF op2 = particles[j]->move(1.0f);
					CParticleF oq2 = particles[j]->next->move(1.0f);
					float va1 = GetVisualAngle2(x.m_X, x.m_Y, oq.m_X, oq.m_Y, op.m_X, op.m_Y);
					float va2 = GetVisualAngle2(op2.m_X, op2.m_Y, oq.m_X, oq.m_Y, op.m_X, op.m_Y);
					float va3 = GetVisualAngle2(x.m_X, x.m_Y, op.m_X, op.m_Y, oq.m_X, oq.m_Y);
					float va4 = GetVisualAngle2(oq2.m_X, oq2.m_Y, op.m_X, op.m_Y, oq.m_X, oq.m_Y);
					bool c1 = cos(va1) > cos(va2); //condition 1
					bool c2 = cos(va3) > cos(va4); //condition 2
					bool c3 = sin(va1)*sin(va2)>0; //condition 3
					bool c4 = sin(va3)*sin(va4)>0; //condition 4
					if(c1 && c2 && c3 && c4)
					{
						t = t0;
						res.first = t;
						res.second = particles[i];
					}
				}
			}
		}
	}
	return res;
}

int
hierarcyConsistencyCheck()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	for(int i=0; i<factory->particles.size(); ++i)
	{
		MovingParticle* p = factory->particles[i];
		if(p->type == Initial)
		{
			vector<MovingParticle*> A;
			vector<MovingParticle*> Q(1, p);
			set<MovingParticle*> S;
			S.insert(p);
			while(Q.empty()==false)
			{
				vector<MovingParticle*> Q2;
				for(int j=0; j<Q.size(); ++j)
				{
					MovingParticle* q = Q[j];
					if(q->type == Axis)
					{
						A.push_back(q);
					}
					for(int k=0; k<q->anscestors.size(); ++k)
					{
						MovingParticle* r = q->anscestors[k];
						if(S.find(r) == S.end())
						{
							Q2.push_back(r);
							S.insert(r);
						}
					}
				}
				Q = Q2;
			}
			printf("%d [%f,%f]->[%f,%f] has %d ancestors.\n", p->id, p->p0.m_X, p->p0.m_Y, p->p.m_X, p->p.m_Y, A.size());
		}
	}
	for(int i=0; i<factory->particles.size(); ++i)
	{
		MovingParticle* p = factory->particles[i];
		if(p->type == Axis)
		{
			vector<MovingParticle*> A;
			vector<MovingParticle*> Q(1, p);
			set<MovingParticle*> S;
			S.insert(p);
			while(Q.empty()==false)
			{
				vector<MovingParticle*> Q2;
				for(int j=0; j<Q.size(); ++j)
				{
					MovingParticle* q = Q[j];
					if(q->type == Initial)
					{
						A.push_back(q);
					}
					for(int k=0; k<q->descendents.size(); ++k)
					{
						MovingParticle* r = q->descendents[k];
						if(S.find(r) == S.end())
						{
							Q2.push_back(r);
							S.insert(r);
						}
					}
				}
				Q = Q2;
			}
			printf("%d [%f,%f]->[%f,%f] has %d descendents.\n", p->id, p->p0.m_X, p->p0.m_Y, p->p.m_X, p->p.m_Y, A.size());
		}
	}
	return 0;
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
quickFinish(vector<MovingParticle*>& particles, float time)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	if(particles.size()==1)
	{
		if(particles[0]->type != Axis)
		{
			MovingParticle* p = factory->makeParticle(particles[0]->p, Axis, time);
			setRelation(particles[0], p);
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
		}
	}
	else if(particles.size()==3)
	{
		pair<float,float> param = _IntersectConvexPolygon::intersect(particles[0]->p0, particles[0]->move(1.0), particles[1]->p0, particles[1]->move(1.0));
		if(param.first != param.first)
		{
			//data are colinear?
		}
		else
		{
			float x1 = (1.0-param.first) * particles[0]->p0.m_X + param.first * particles[0]->move(1.0).m_X;
			float y1 = (1.0-param.first) * particles[0]->p0.m_Y + param.first * particles[0]->move(1.0).m_Y;
			CParticleF f(x1,y1);
			MovingParticle* pnew = factory->makeParticle(f, Axis, time);
			for(int i=0; i<3; ++i)
			{
				particles[i]->p=f;
				setRelation(particles[i], pnew);
			}
		}
	}
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

/*
remove duplicates and kinks. 
*/
void 
removeDegenerate(vector<MovingParticle*>& particles, float time, float eps)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<MovingParticle*> setv;
	while(particles.size()>1) //at least two points are needed. else particles can vanish.
	{
		bool bchanged = false;
		for(int i=0; i<particles.size(); ++i)
		{
			MovingParticle* p = particles[i];
			MovingParticle* q = p->next;
			MovingParticle* r = p->prev;
			float dpq = Distance(p->p, q->p);
			float dpr = Distance(p->p, r->p);
			float dqr = Distance(q->p, r->p);
			if(dpq< eps) //points are close enough
			{
				//remove the two identical particles and create a new one at the same position.
				MovingParticle* pnew = factory->makeParticle(p->p, Collide, time);
				particles[i] = pnew;
				setRelation(p, particles[i]);
				setRelation(q, particles[i]);
				particles[i]->next = q->next;
				particles[i]->prev = p->prev;
				particles[i]->next->prev = particles[i];
				particles[i]->prev->next = particles[i];
				int i2 = (i + 1) % particles.size();
				particles.erase(particles.begin()+i2);
				bchanged = true;
				setv.push_back(pnew); // velocity needs to be set
				break;
			}
			else if(particles.size()>=3)
			{
				pair<float,float> param = _IntersectConvexPolygon::intersect(p->p0, p->p, q->p0, q->p);
				//just in case when particles cross over between p and q.
				if(param.first>0 && param.first <=1.0 && param.second>0 && param.second<=1.0)
				{
					float x = (1.0-param.first)*p->p0.m_X + param.first*p->p.m_X;
					float y = (1.0-param.first)*p->p0.m_Y + param.first*p->p.m_Y;
					particles[i] = factory->makeParticle(CParticleF(x, y), Collide, time);
					setv.push_back(particles[i]); // velocity needs to be set
					setRelation(p, particles[i]);
					setRelation(q, particles[i]);
					particles[i]->next = q->next;
					particles[i]->prev = p->prev;
					particles[i]->next->prev = particles[i];
					particles[i]->prev->next = particles[i];
					int i2 = (i + 1) % particles.size();
					particles.erase(particles.begin()+i2);
					bchanged = true;
					break;
				}
				//two neighbors are too close -> angle too small
				//else if(dpr > eps && GetVisualAngle(q->p.m_X, q->p.m_Y, r->p.m_X, r->p.m_Y, p->p.m_X, p->p.m_Y) < 0.001)
				else if(dpr > eps && (Distance2LineSegment(p->p, q->p, r->p)<eps || Distance2LineSegment(p->p, r->p, q->p)<eps))
				{
					MovingParticle* z = factory->makeParticle(p->p, Collide, time);
					//MovingParticle* z = factory->makeParticle(p->p, Merge);
					//MovingParticle* z = factory->makeParticle(p->p, Axis);
					if(dpq < dpr)//move it to the closest neighbor
					{
						z->p = q->p;
						setRelation(p, z);
						setRelation(q, z);
					}
					else
					{
						z->p = r->p;
						setRelation(p, z);
						setRelation(r, z);
					}
					//one of the neighbor (the closest one) is at the same location.
					//thus, they will be replace by a new particle on the next iteration after 'break'
					z->next = q;
					z->prev = r;
					z->next->prev = z;
					z->prev->next = z;
					particles[i] = z;
					bchanged = true;
					break;
				}
			}
		}
		if(bchanged == false) break;
	}
	//set velocity of new particles after all degenerates are removed.
	for(int i=0; i<setv.size(); ++i)
	{
		setv[i]->setVelocity();
	}
}

/*
Go through a chanin of points (in PARTICLES) with possible splits, separate it into a set of polygons.
Note that polygons can include a degenerate one with only two points whena a split spans multiple points (where a skeleton
is formed.)
CMAP gives pairs of points at splits. If cmap[x]==y, then cmap[y]==x. If it is not a split point, then cmap[x] is NULL.
The procedure works in the following manner.
1. for each split point (denoted by non-null CMAP), create a new particle.
2. if possible, complete a polygon, and add the polygon to a vector.
3. else, add the partial trace to the stack.
4. after going through all particles, stitch partial traces in the stack to complete the last polygon.
5. for each NON-MERGE polygon, set neighbors/ancestors/descendents/velocity of new particles (indicated by NULL neighbors).
6. for each MERGE (polygon - i.e. a group comprised of only two points), set ancestors/descendents. We need to add the twin particle and its 
partner to its ascendent list.
*/
vector<vector<MovingParticle*>>
tracePolygon(vector<MovingParticle*>& particles,
	map<MovingParticle*,MovingParticle*>& cmap,
	float time)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<vector<MovingParticle*>> polygons;
	stack<vector<MovingParticle*>> partials; //store partially formed polygons
	vector<MovingParticle*> trace;
	map<MovingParticle*,pair<MovingParticle*,MovingParticle*>> children; //given the new particle, provides a pair of splits
	map<MovingParticle*,MovingParticle*> addition; //newly added particle given the original.
	for(int i=0; i<particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		if(cmap[p] == NULL)
		{
			trace.push_back(p);
		}
		else
		{
			MovingParticle* q = factory->makeParticle(p->p, Unknown, time);
			children[q] = pair<MovingParticle*,MovingParticle*>(p, cmap[p]); //for setting anscestory later on
			addition[p] = q; //used to associate a nearby non-Merge split. 
			if(trace.empty())
			{
			}
			else if(cmap[p] == children[trace[0]].first) //complete the polygon
			{
				polygons.push_back(trace);
			}
			else if(partials.empty()==false && cmap[p] == children[partials.top()[0]].first) //complete with the stack
			{
				vector<MovingParticle*> top = partials.top();
				top.insert(top.end(), trace.begin(), trace.end());
				polygons.push_back(top);
				partials.pop(); //remove the top 
			}
			else //does not complete.
			{
				partials.push(trace);
			}
			trace.clear();
			trace.push_back(q);
		}
	}
	//the stack contains a remaining polygon.
	while(partials.empty() == false)
	{
		vector<MovingParticle*> tr = partials.top();
		trace.insert(trace.begin(), tr.begin(), tr.end());
		partials.pop();
	}
	if(trace.empty()==false)
	{
		polygons.push_back(trace);
	}
	for(int i=0; i<polygons.size(); ++i)
	{
		int n = polygons[i].size();
		for(int j=0; j<n; ++j)
		{
			MovingParticle* p = polygons[i][j];
			if(p->type == Unknown) //new particle
			{
				MovingParticleType type = n==2 ? Merge: Split;
				p->type = type;
				int j0 = (j-1+n) % n;
				int j2 = (j+1) % n;
				setNeighbors(p, polygons[i][j0], polygons[i][j2]);
				//if(children[p].first->type != Dummy)
				{
					setRelation(children[p].first, p);
				}
				//if(children[p].second->type != Dummy)
				{
					setRelation(children[p].second, p);
				}
				if(type == Split)
				{
					p->setVelocity();
				}
				if(type == Merge)
				{
					//The next line creates bidirectional relation between two end points.
					//This was to prevent any gaps from froming when polygons are split per axis point.
					//Now, I an investigating smoothing and subsequent grouping, which is drastically different from
					//the splitting approach. Thus, disabled this line.
					//setRelation(p, polygons[i][j2]);

					//find the origin of its counter-part
					//I don't remember what this code is for... Do I still need it?
					//MovingParticle* p0 = children[p].second;
					//setRelation(p, addition[p0]);
				}
			}
		}
	}
	return polygons;
}

/*
Loop through the current polygon, find any points that have overlaps (either splitting by a concave edge or merging by parallel edges),
and collect entry and exiting points of each overlap. 
For each entry and exiting point, provide a bridging particle to be inserted.
After the bridges are inserted, trace the particles once more and divide it into multiple polygons.
*/
vector<DeformingPolygon>
fixPolygon(DeformingPolygon& polygon, float time, float eps, bool debugOn)
{
	if(debugOn)
	{
		printf("Before removing degenerate.\n");
		for(int i=0; i<polygon.particles.size(); ++i)
		{
			MovingParticle* p = polygon.particles[i];
			printf("%d %f %f %f %f %f %f %d %d\n", i+1, p->p0.m_X, p->p0.m_Y, p->p.m_X, p->p.m_Y, p->v[0], p->v[1], p->id, p->type);
		}
	}
	removeDegenerate(polygon.particles, time, eps);
	if(debugOn)
	{
		printf("After removing degenerate.\n");
		for(int i=0; i<polygon.particles.size(); ++i)
		{
			MovingParticle* p = polygon.particles[i];
			printf("%d %f %f %f %f %f %f %d %d\n", i+1, p->p0.m_X, p->p0.m_Y, p->p.m_X, p->p.m_Y, p->v[0], p->v[1], p->id, p->type);
		}
	}
	if(polygon.particles.size() <= 3)
	{
		time += 0;
		//return vector<DeformingPolygon>();
	}

	ParticleFactory* factory = ParticleFactory::getInstance();
	int N = polygon.particles.size();
	map<MovingParticle*,MovingParticle*> cmap;
	vector<MovingParticle*> dummy;

	//First check if any overlap and place dummy if necessary.
	for(int i=0; i<N; ++i)
	{
		cmap[polygon.particles[i]] = NULL;
	}
	for(int i=0; i<N; ++i)
	{
		MovingParticle* p = polygon.particles[i];
		for(int j=0; j<N; ++j)
		{
			MovingParticle* q = polygon.particles[j];
			if(p==q || q->next==p) continue;
			float sep = Distance(p->p,q->p);
			float sep2 = Distance(p->p,q->next->p);
			if(sep < eps) 
			{
				cmap[p] = q;
				cmap[q] = p;
			}
			else if(sep2 < eps)
			{
				cmap[p] = q->next;
				cmap[q->next] = p;
			}
			else
			{
				CParticleF y;
				if(overlap(p->p, q->p, q->next->p, p->p0, y, eps))//needs a dummy...
				{
					MovingParticle* x = factory->makeParticle(y, Dummy, time);
					//setRelation(x0, x);
					setNeighbors(x, q, q->next);
					x->setVelocity();
					dummy.push_back(x);
					cmap[p] = x;
					cmap[x] = p;
				}
			}
			if(cmap[p] != NULL) break;
		}
	}

	vector<MovingParticle*> particles = vectorize(polygon.particles[0]);

	if(debugOn)
	{
		printf("After insertion of dummies.\n");
		for(int i=0; i<particles.size(); ++i)
		{
			MovingParticle* p = particles[i];
			printf("%d %f %f %f %f %f %f %d %d\n", i+1, p->p0.m_X, p->p0.m_Y, p->p.m_X, p->p.m_Y, p->v[0], p->v[1], p->id, p->type);
		}
	}
	vector<vector<MovingParticle*>> polygons = tracePolygon(particles, cmap, time);
	if(debugOn)
	{
		printf("After tracing and splitting polygons.\n");
		for(int j=0; j<polygons.size(); ++j)
		{
			printf("Polygon %d\n", j+1);
			for(int i=0; i<polygons[j].size(); ++i)
			{
				MovingParticle* p = polygons[j][i];
				printf("%d %f %f %f %f %f %f %d %d\n", i+1, p->p0.m_X, p->p0.m_Y, p->p.m_X, p->p.m_Y, p->v[0], p->v[1], p->id, p->type);
			}
		}
	}

	vector<DeformingPolygon> shapes;
	if(polygons.size() == 1)
	{
		polygon.particles = polygons[0];
		shapes.push_back(polygon);
	}
	else
	{
		for(int i=0; i<polygons.size(); ++i)
		{
			DeformingPolygon dp(time, polygon.id);
			dp.particles = polygons[i];
			shapes.push_back(dp);
		}
	}
	return shapes;
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
		poly[i]->p = poly[i]->move(0.1);
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


mxArray*
StoreStraightAxes(vector<DeformingPolygon>& polygons) 
{
	set<MovingParticle*> pset; //collect a set of particles in the representation
	for(int i=0; i<polygons.size(); ++i)
	{
		for(int j=0; j<polygons[i].particles.size(); ++j)
		{
			pset.insert(polygons[i].particles[j]);
			for(int k=0; k<polygons[i].particles[j]->descendents.size(); ++k)
			{
				pset.insert(polygons[i].particles[j]->descendents[k]);
			}
			for(int k=0; k<polygons[i].particles[j]->anscestors.size(); ++k)
			{
				pset.insert(polygons[i].particles[j]->anscestors[k]);
			}
		}
	}
	//sort unique set of particles based on their id numbers
	vector<pair<float,MovingParticle*>> pairs;
	for(set<MovingParticle*>::iterator it=pset.begin(); it != pset.end(); it++)
	{
		MovingParticle* p = *it;
		pairs.push_back(pair<float,MovingParticle*>(p->id, p));
	}
	sort(pairs.begin(), pairs.end());

	//map particles to their indices in the vector 
	map<MovingParticle*,int> imap;
	for(int i=0; i<pairs.size(); ++i)
	{
		imap[pairs[i].second] = i;
	}

	mxArray* plhs;
	const mwSize dimsS[2]={1, 1};
	const char* fieldNames[] = {"particles", "polygons"};
	const int NumberOfFields = sizeof(fieldNames)/sizeof(fieldNames[0]);
	/* Create a structs. */ 
	plhs = mxCreateStructArray(2, dimsS, NumberOfFields, fieldNames);

	// Store particles
	{
		// Create a struct array for particles
		int n = imap.size();
		const mwSize dims[2]={n, 1};
		const char* fieldNames2[] = {"p0", "p", "v", "type", "id", "created", "prev", "next", "ancestors", "descendents"};
		const int NumberOfFields2 = sizeof(fieldNames2)/sizeof(fieldNames2[0]);
		mxArray* mxPtrParticles = mxCreateStructArray(2, dims, NumberOfFields2, fieldNames2);
		for(int i=0; i<pairs.size(); ++i)
		{
			MovingParticle* p = pairs[i].second;
			{
				//p0
				const mwSize dims2[2] = {1, 2};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxSINGLE_CLASS, mxREAL);
				float* fptr = (float*) mxGetData(mxptr);
				fptr[0] = p->p0.m_X;
				fptr[1] = p->p0.m_Y;
				mxSetField(mxPtrParticles, i, fieldNames2[0], mxptr);
			}
			{
				//p
				const mwSize dims2[2] = {1, 2};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxSINGLE_CLASS, mxREAL);
				float* fptr = (float*) mxGetData(mxptr);
				fptr[0] = p->p.m_X;
				fptr[1] = p->p.m_Y;
				mxSetField(mxPtrParticles, i, fieldNames2[1], mxptr);
			}
			{
				//v
				const mwSize dims2[2] = {1, 2};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxSINGLE_CLASS, mxREAL);
				float* fptr = (float*) mxGetData(mxptr);
				fptr[0] = p->v[0];
				fptr[1] = p->v[1];
				mxSetField(mxPtrParticles, i, fieldNames2[2], mxptr);
			}
			{
				//type
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				fptr[0] = p->type;
				mxSetField(mxPtrParticles, i, fieldNames2[3], mxptr);
			}
			{
				//id
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				fptr[0] = p->id;
				mxSetField(mxPtrParticles, i, fieldNames2[4], mxptr);
			}
			{
				//creation time
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxSINGLE_CLASS, mxREAL);
				float* fptr = (float*) mxGetData(mxptr);
				fptr[0] = p->created;
				mxSetField(mxPtrParticles, i, fieldNames2[5], mxptr);
			}
			{
				//prev
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				if(p->prev == NULL)
				{
					fptr[0] = 0;
				}
				else
				{
					fptr[0] = imap[p->prev] + 1;
				}
				mxSetField(mxPtrParticles, i, fieldNames2[6], mxptr);
			}
			{
				//next
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				if(p->next == NULL)
				{
					fptr[0] = 0;
				}
				else
				{
					fptr[0] = imap[p->next] + 1;
				}
				mxSetField(mxPtrParticles, i, fieldNames2[7], mxptr);
			}
			{
				//anscestors
				const mwSize dims2[2] = {1, p->anscestors.size()};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				for(int j=0; j<p->anscestors.size(); ++j)
				{
					fptr[j] = imap[p->anscestors[j]] + 1;
				}
				mxSetField(mxPtrParticles, i, fieldNames2[8], mxptr);
			}
			{
				//descendents
				const mwSize dims2[2] = {1, p->descendents.size()};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				for(int j=0; j<p->descendents.size(); ++j)
				{
					fptr[j] = imap[p->descendents[j]] + 1;
				}
				mxSetField(mxPtrParticles, i, fieldNames2[9], mxptr);
			}
		}
		mxSetField(plhs, 0, fieldNames[0], mxPtrParticles);
	}
	//store polygons
	{
		int n = polygons.size();
		const mwSize dims[2]={n, 1};
		const char* fieldNames2[] = {"created", "id", "parent_id", "particles"};
		const int NumberOfFields2 = sizeof(fieldNames2)/sizeof(fieldNames2[0]);
		mxArray* mxPtrParticles = mxCreateStructArray(2, dims, NumberOfFields2, fieldNames2);
		for(int i=0; i<n; i++)
		{
			{
				//creation time
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxSINGLE_CLASS, mxREAL);
				float* fptr = (float*) mxGetData(mxptr);
				fptr[0] = polygons[i].created;
				mxSetField(mxPtrParticles, i, fieldNames2[0], mxptr);
			}
			{
				//id
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				fptr[0] = polygons[i].id;
				mxSetField(mxPtrParticles, i, fieldNames2[1], mxptr);
			}
			{
				//parent_id
				const mwSize dims2[2] = {1, 1};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				fptr[0] = polygons[i].parent_id;
				mxSetField(mxPtrParticles, i, fieldNames2[2], mxptr);
			}
			{
				//particles
				int m = polygons[i].particles.size();
				const mwSize dims2[2] = {1, m};
				mxArray* mxptr = mxCreateNumericArray(2, dims2, mxINT32_CLASS, mxREAL);
				int* fptr = (int*) mxGetData(mxptr);
				for(int j=0; j<m; ++j)
				{
					fptr[j] = imap[polygons[i].particles[j]] + 1; //one based
				}
				mxSetField(mxPtrParticles, i, fieldNames2[3], mxptr);
			}
		}
		mxSetField(plhs, 0, fieldNames[1], mxPtrParticles);
	}
	return plhs;
}

bool LoadStraightAxes(vector<DeformingPolygon>& polygons, const mxArray *prhs)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	//get particles
	mxArray* pa = mxGetField(prhs, 0, "particles");
	mxArray* pb = mxGetField(prhs, 0, "polygons");
	vector<MovingParticle*> particles;
	{
		int n = mxGetNumberOfElements(pa);
		for(int i=0; i<n; ++i)
		{
			CParticleF p0, p;
			float v[2];
			MovingParticleType type;
			int id;
			float time;
			{
				mxArray* paa = mxGetField(pa, i, "p0");
				float* fptr = (float*) mxGetData(paa);
				p0.m_X  = fptr[0];
				p0.m_Y = fptr[1];
			}
			{
				mxArray* paa = mxGetField(pa, i, "p");
				float* fptr = (float*) mxGetData(paa);
				p.m_X  = fptr[0];
				p.m_Y = fptr[1];
			}
			{
				mxArray* paa = mxGetField(pa, i, "v");
				float* fptr = (float*) mxGetData(paa);
				v[0]  = fptr[0];
				v[1] = fptr[1];
			}
			{
				mxArray* paa = mxGetField(pa, i, "type");
				int* fptr = (int*) mxGetData(paa);
				int t = (int) fptr[0];
				switch(t)
				{
				case 0:
					type = Unknown;
					break;
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
			}
			{
				mxArray* paa = mxGetField(pa, i, "id");
				int* fptr = (int*) mxGetData(paa);
				id  = (int) fptr[0];
			}
			{
				mxArray* paa = mxGetField(pa, i, "created");
				float* fptr = (float*) mxGetData(paa);
				time  = (float) fptr[0];
			}
			MovingParticle* particle = factory->makeParticle(p0, type, time);
			particle->p = p;
			//particle->id = id;
			particle->v[0] = v[0];
			particle->v[1] = v[1];
			particles.push_back(particle);
		}
		//now go through again and set neighbors and ancestors/descendents
		for(int i=0; i<n; ++i)
		{
			{
				mxArray* paa = mxGetField(pa, i, "prev");
				int* fptr = (int*) mxGetData(paa);
				int j  = (int) fptr[0] - 1; //0 based
				if(j>=0)
				{
					particles[i]->prev = particles[j];
				}
				else
				{
					particles[i]->prev = NULL;
				}
			}
			{
				mxArray* paa = mxGetField(pa, i, "next");
				int* fptr = (int*) mxGetData(paa);
				int j  = (int) fptr[0] - 1; //0 based
				if(j>=0)
				{
					particles[i]->next = particles[j];
				}
				else
				{
					particles[i]->next = NULL;
				}
			}
			{
				mxArray* paa = mxGetField(pa, i, "ancestors");
				int* fptr = (int*) mxGetData(paa);
				int m = mxGetNumberOfElements(paa);
				for(int k=0; k<m; ++k)
				{
					int j  = (int) fptr[k] - 1; //0 based
					particles[i]->anscestors.push_back(particles[j]);
				}
			}
			{
				mxArray* paa = mxGetField(pa, i, "descendents");
				int* fptr = (int*) mxGetData(paa);
				int m = mxGetNumberOfElements(paa);
				for(int k=0; k<m; ++k)
				{
					int j  = (int) fptr[k] - 1; //0 based
					particles[i]->descendents.push_back(particles[j]);
				}
			}
		}
	}
	//poygons
	{
		int n = mxGetNumberOfElements(pb);
		for(int i=0; i<n; ++i)
		{
			float time;
			{
				mxArray* paa = mxGetField(pb, i, "created");
				float* fptr = (float*) mxGetData(paa);
				time  = fptr[0];
			}
			int id;
			{
				mxArray* paa = mxGetField(pb, i, "id");
				int* fptr = (int*) mxGetData(paa);
				id  = (int)fptr[0];
			}
			int parent_id;
			{
				mxArray* paa = mxGetField(pb, i, "parent_id");
				int* fptr = (int*) mxGetData(paa);
				parent_id  = (int)fptr[0];
			}
			vector<MovingParticle*> ps;
			vector<int> vidx;
			{
				mxArray* paa = mxGetField(pb, i, "particles");
				int m = mxGetNumberOfElements(paa);
				int* fptr = (int*) mxGetData(paa);
				for(int j=0; j<m; ++j)
				{
					ps.push_back(particles[fptr[j]-1]); //0 index
					vidx.push_back(fptr[j]);
				}
			}
			DeformingPolygon poly(time);
			poly.id = id;
			poly.parent_id = parent_id;
			poly.particles = ps;
			polygons.push_back(poly);
		}
	}
	return true;
}

//	int MovingParticle::_id = 0;
//	int DeformingPolygon::_id = 0;

}