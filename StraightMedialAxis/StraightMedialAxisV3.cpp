#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
#include <set>
#include <hash_map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <DisjointSet.h>
#include <TriangulationWithMemento.h>
#include <FragmentInfo.h>
//#include <TriangulationHelper.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
typedef CParticleF graphKey;

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;

template<class X, class Y, class Z>
struct triple
{
	triple(X x, Y y, Z z)
	{
		first = x;
		second = y;
		third = z;
	}
	X first;
	Y second;
	Z third;
};

/*
This method finds an intersection point between a plane and a line. The plane is defined by op and z where the former is a point on the plane
and the latter is a vector perpendicular to the plane. The line is defined by ol and u where the former is a point on the line and
the latter is a vector along the line.
It returns a parameter t, such that the intersection point is at ol + t * u.
*/
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
	vx = vx / len0;
	vy = vy / len0;
	return CParticleF(vx, vy);
}

/*
Check for each point, if it is convex or not with respect to the polygon defined by points[].
*/
vector<bool>
convexity(const vector<CParticleF>& points)
{
	vector<bool> bflag(points.size());
	int cw = ClockWise(points);
	for(int i=0; i<points.size(); ++i)
	{
		int j=(i-1+points.size()) % points.size();
		int k=(i+1) % points.size();
		int cw2 = ClockWise(points[j].m_X, points[j].m_Y, points[i].m_X, points[i].m_Y, points[k].m_X, points[k].m_Y);
		
		if(cw * cw2 < 0)
		{
			bflag[i] = false;
		}
		else
		{
			bflag[i] = true;
		}
	}
	return bflag;
}

struct MovingParticle
{
	MovingParticle(CParticleF& o, CParticleF& a, CParticleF& b, bool convex, int id=0)
	{
		this->p = o;
		this->convex = convex;
		setVelocity(o, a, b, convex);
		this->id = id;
	}
	void setVelocity(CParticleF& o, CParticleF& a, CParticleF& b, bool convex)
	{
		CParticleF bs = bisector(o, a, b);

		float ang = (PI - GetVisualAngle(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y)) / 2.0;
		float cs = cos(ang);
		/*float len = 1.0f;
		if(Abs(cs) > .05)
		{
			len = 1.0 / cs;
		}*/
		float len = 1.0f / cos(ang);
		bs.m_X *= len;
		bs.m_Y *= len;
		if(convex == false)
		{
			bs.m_X = -bs.m_X;
			bs.m_Y = -bs.m_Y;
		}
		this->v = bs;
	}
	CParticleF move(float t)
	{
		p.m_X = p.m_X + t * v.m_X;
		p.m_Y = p.m_Y + t * v.m_Y;
		return p;
	}

	CParticleF p;
	CParticleF v;
	bool convex;
	int id;
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
	MovingParticle* makeParticle(CParticleF& p, CParticleF& prev, CParticleF& next, bool convex)
	{
		MovingParticle* particle = new MovingParticle(p, prev, next, convex, _id);
		_id++;
		particles.push_back(particle);
		return particle;
	}
private:
	static ParticleFactory* _instance;
	static int _id;
	ParticleFactory(){};
	vector<MovingParticle*> particles;
};

ParticleFactory* ParticleFactory::_instance = NULL;
int ParticleFactory::_id = 0;

vector<MovingParticle*>
createMovingParticles(vector<CParticleF>& points)
{
	vector<MovingParticle*> particles;
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<bool> bc = convexity(points);
	for(int i=0; i<points.size(); ++i)
	{
		int j=(i-1+points.size()) % points.size();
		int k=(i+1) % points.size();
		MovingParticle* p = factory->makeParticle(points[i], points[j], points[k], bc[i]);
		particles.push_back(p);
	}
	return particles;
}

struct EventDescription
{
	float t;
	MovingParticle* p; //event particle.
	MovingParticle* q; //a colliding particle if edge event. the first particle of the colliding segment if split event.
	CParticleF loc; //the location of the event
	//bool edgeEvent;
	bool operator <(EventDescription ev)
	{
		return t < ev.t;
	}
};

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

vector<CParticleF>
removeDegenerate(vector<CParticleF> points)
{
	while(true)
	{
		vector<CParticleF> result;
		bool bChanged = false;
		for(int i=0; i<points.size(); ++i)
		{
			int i0=(i-1+points.size()) % points.size();
			int i2=(i+1) % points.size();
			if(coLinear(points[i0], points[i], points[i2]))
			{
				bChanged = true;
			}
			else
			{
				result.push_back(points[i]);
			}
		}
		if(bChanged == false)
		{
			return result;
		}
		else
		{
			points = result;
		}
	}
}

vector<vector<MovingParticle*>> 
restructureShape(vector<MovingParticle*>& particles)
{
	vector<vector<MovingParticle*>> splitted;
	//1. remove duplicate points
	float eps = 0.1;
	while(true)
	{
		bool bchanged = false;
		for(int i=0; i<particles.size()-1; ++i)
		{
			int j=i+1;
			if(Abs(particles[i]->p.m_X-particles[j]->p.m_X)<eps && Abs(particles[i]->p.m_Y-particles[j]->p.m_Y)<eps)
			{
				particles.erase(particles.begin()+j);
				bchanged = true;
				break;
			}
		}
		if(bchanged==false) break;
	}
	//2. detect colinear segments
	{
		//first cluster colinear segments together
		vector<Node<int>*> nodes;
		for(int i=0; i<particles.size(); ++i)
		{
			nodes.push_back(makeset(i));
		}
		for(int i=0; i<particles.size(); ++i)
		{
			int i2 = (i+1) % particles.size();
			for(int j=i+1; j<particles.size(); ++j)
			{
				int j2 = (j+1) % particles.size(); 
				if(coLinear(particles[i]->p, particles[i2]->p, particles[j]->p, particles[j2]->p))
				{
					merge(nodes[i], nodes[j]);
				}
			}
		}
		vector<Node<int>*> cl = clusters(nodes);
		//figure out which segments are colinear.
		map<Node<int>*,int> countmap;
		for(int i=0; i<cl.size(); ++i) 
		{
			countmap[cl[i]] = 0;
		}
		for(int i=0; i<nodes.size(); ++i) 
		{
			countmap[findset(nodes[i])]++;
		}
		vector<bool> bkeep(nodes.size(), false);
		int start = -1;
		for(int i=0; i<nodes.size(); ++i)
		{
			if(countmap[findset(nodes[i])]==1)
			{
				bkeep[i] = true;
			}
			else
			{
				start = i;
			}
		}
		if(start < 0) //there was no colinear segment. No splitting is required.
		{
			splitted.push_back(particles);
		}
		else
		{
			vector<MovingParticle*> shape;
			for(int i=0; i<particles.size(); ++i)
			{
				int j = (start + 1 + i) % particles.size(); 
				if(bkeep[i])
				{
					shape.push_back(particles[j]);
				}
				else
				{
					if(shape.size() > 2)
					{
						splitted.push_back(shape);
						shape.clear();
					}
				}
			}
		}
	}
	//3. find edge-splitting vertices (split event) and insert a new vertex.
	vector<vector<MovingParticle*>> splitted2;
	for(int s=0; s<splitted.size(); ++s)
	{
		vector<MovingParticle*> shape = splitted[s];
		vector<pair<MovingParticle*,MovingParticle*>> pairs;
		for(int i=0; i<shape.size(); ++i)
		{
			for(int j=0; j<shape.size(); ++j)
			{
				int j0 = (j-1+shape.size()) % shape.size();
				int j2 = (j+1) % shape.size();
				if(j==i || j==j0 || j==j2) continue;
				if(coLinear(shape[j]->p, shape[j2]->p, shape[i]->p))
				{
					pair<MovingParticle*,MovingParticle*> p(shape[i], shape[j]);
					pairs.push_back(p);
					break;
				}
			}
		}
		if(pairs.empty()) //no vertex splitted the shape.
		{
			splitted2.push_back(shape);
		}
		else
		{
			ParticleFactory* factory = ParticleFactory::getInstance();
			vector<pair<MovingParticle*,MovingParticle*>> clones;
			for(int i=0; i<pairs.size(); ++i)
			{
				int k = distance(shape.begin(), find(shape.begin(), shape.end(), pairs[i].first));
				MovingParticle* prev = shape[(k-1+shape.size()) % shape.size()];
				MovingParticle* next = shape[(k-1+shape.size()) % shape.size()];
				MovingParticle* pc = factory->makeParticle(pairs[i].first->p, prev->p, next->p, true);
				int k2 = distance(shape.begin(), find(shape.begin(), shape.end(), pairs[i].second));
				k2 = (k2 + 1) % shape.size();
				shape.insert(shape.begin()+k2, pc);
				pair<MovingParticle*,MovingParticle*> clone(pairs[i].first, pc);
				clones.push_back(clone);
			}
			for(int i=0; i<clones.size(); ++i)
			{
				vector<MovingParticle*> loop;
				MovingParticle* p = clones[i].first;
				int k = distance(shape.begin(), find(shape.begin(), shape.end(), p));
				while(p != clones[i].second)
				{
					loop.push_back(p);
					k = (k + 1) % shape.size(); 
					p = shape[k];
				}
			}
		}
	}
}


vector<EventDescription>
futureEdgeEvents(vector<MovingParticle*>& v)
{
	vector<EventDescription> evs;
	for(int i=0; i<v.size(); ++i)
	{
		//if(v[i]->convex == false) continue;
		CParticleF qi(v[i]->p.m_X+v[i]->v.m_X, v[i]->p.m_Y+v[i]->v.m_Y);
		//for(int j=i+1; j<v.size(); ++j)
		int j = (i+1) % v.size();
		{
			//if(v[j]->convex == false) continue;
			if(i+1==10 && j+1==11)
				i+=0;
			CParticleF qj(v[j]->p.m_X+v[j]->v.m_X, v[j]->p.m_Y+v[j]->v.m_Y);
			pair<float,float> param = _IntersectConvexPolygon::intersect(v[i]->p, qi, v[j]->p, qj);
			float t = Max(param.first, param.second);
			if(t>=0)
			{
				EventDescription ev;
				ev.t = t;
				ev.p = v[i];
				ev.q = v[j];
				ev.loc = CParticleF((1-param.first)*v[i]->p.m_X+param.first*qi.m_X, (1-param.first)*v[i]->p.m_Y+param.first*qi.m_Y); 
				evs.push_back(ev);
			}
		}
	}
	sort(evs.begin(), evs.end()); //sort events based on their time stamps
	return evs;
}

/*
For the event to occur, a non-convex particle touches moving segment that is not adjacent to the particle.
To detect such event, we consider the moving line segment as a plane with the third axis as time.
We then compute the time the particle (in 3D) intersects the plane.
We then need to check if the time is valid. It has to be forward looking from both particle and intersecting line segment.
It has to be within a sweep of the line segment.
*/
EventDescription
nextSplitEvent(vector<MovingParticle*>& v)
{
	EventDescription ev;
	ev.t = std::numeric_limits<float>::infinity();
	for(int i=0; i<v.size(); ++i)
	{
		if(v[i]->convex) continue;
		CParticleF ol = v[i]->p;
		CParticleF ul = v[i]->v;
		ul.m_Z = 1.0;
		for(int j=0; j<v.size(); ++j)
		{
			int j0 = (j-1+v.size()) % v.size();
			int j2 = (j+1) % v.size();
			if(i==j || i==j2 || i==j0) continue;
			//if(i==j) continue;
			CParticleF op = v[j]->p;
			CParticleF oq = v[j2]->p;
			CParticleF ov = NormalizedDirection(op, oq);
			CParticleF up0(ov.m_Y, -ov.m_X, 1.0f);
			CParticleF up = crossProduct(ov, up0);
			if(up.m_X*v[j]->v.m_X + up.m_Y*v[j]->v.m_Y < 0)
			{
				up.m_X = -up.m_X;
				up.m_Y = -up.m_Y;
			}
			float t = intersectPlaneAndLine(op, up, ol, ul);
			CParticleF loc(ol.m_X+t*ul.m_X, ol.m_Y+t*ul.m_Y, ol.m_Z+t*ul.m_Z);
			if(t >= 0 && t<ev.t) 
			{
				//check if the intersection point is within the sweep of the line segment.
				CParticleF op2(op.m_X + t * v[j]->v.m_X, op.m_Y + t * v[j]->v.m_Y);
				CParticleF oq2(oq.m_X + t * v[j2]->v.m_X, oq.m_Y + t * v[j2]->v.m_Y);
				if(Distance2LineSegment(op2,  oq2, loc)<1.0f)
				{
					ev.t = t;
					ev.p = v[i];
					ev.q = v[j];
					ev.loc = CParticleF(ol.m_X + t*ul.m_X, ol.m_Y + t*ul.m_Y); 
				}
			}
		}
	}
	return ev;
}

pair<vector<MovingParticle*>,vector<MovingParticle*>>
splitPolygon(vector<MovingParticle*>& v, 
			 EventDescription& ev)
{
	int k1 = distance(v.begin(), find(v.begin(), v.end(), ev.q));
	int k0 = (k1-1+v.size()) % v.size();
	int k2 = (k1+1) % v.size();
	ParticleFactory* factory = ParticleFactory::getInstance();
	MovingParticle* z = factory->makeParticle(v[k1]->p, v[k0]->p, v[k2]->p, v[k1]->convex);
	v.insert(v.begin()+k2, z);
	vector<MovingParticle*> poly1;
	{
		int k = distance(v.begin(), find(v.begin(), v.end(), ev.p));
		MovingParticle* p = v[k];
		while(p != z)
		{
			poly1.push_back(p);
			k = (k + 1) % v.size(); 
			p = v[k];
		}
	}
	vector<MovingParticle*> poly2;
	{
		int k = distance(v.begin(), find(v.begin(), v.end(), z));
		MovingParticle* p = v[k];
		while(p != ev.p)
		{
			poly2.push_back(p);
			k = (k + 1) % v.size(); 
			p = v[k];
		}
	}
	pair<vector<MovingParticle*>,vector<MovingParticle*>> splitted;
	splitted.first = poly1;
	splitted.second = poly2;
	return splitted;
}

void 
StraightAxisPropagation(vector<MovingParticle*> particles,
						map<MovingParticle*,Vertex<CParticleF>*>& hmap)
{
	GraphFactory<CParticleF>* factory = GraphFactory<CParticleF>::GetInstance();
	if(particles.size() > 2)
	{
		vector<EventDescription> edgeev = futureEdgeEvents(particles);
		EventDescription splitev = nextSplitEvent(particles);
		float t = splitev.t;
		bool bedge = false;
		if(edgeev.empty() == false && edgeev[0].t < splitev.t)
		{
			bedge = true;
			t = edgeev[0].t;
		}
		if(t >= std::numeric_limits<float>::infinity()) 
		{
			printf("Premature ending...\n");
			for(int j=0; j<particles.size(); ++j)
			{
				printf("%d: (%3.3f,%3.3f), [%3.3f,%3.3f]\n", 
					j+1, particles[j]->p.m_X, particles[j]->p.m_Y, particles[j]->v.m_X, particles[j]->v.m_Y);
			}
			return;
		}
		//move the particles by the earliest nest event time
		for(int j=0; j<particles.size(); ++j)
		{
			Vertex<CParticleF>* vt = hmap[particles[j]];
			particles[j]->move(t);
			Vertex<CParticleF>* vt2 = factory->makeVertex(particles[j]->p);
			factory->makeEdge(vt, vt2, 1.0f);
			hmap[particles[j]] = vt2;
		}
		//connect new polygon with edges
		for(int j=0; j<particles.size(); ++j)
		{
			int j2 = (j+1) % particles.size();
			factory->makeEdge(hmap[particles[j]], hmap[particles[j2]], bedge ? 0.0f: 2.0f);
		}
		
		vector<vector<MovingParticle*>> shapes = restructureShape(particles);
		for(int i=0; i<shapes.size(); ++i)
		{
			StraightAxisPropagation(shapes[i], hmap);
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "StraightMedialAxis", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = StraightMedialAxis(P, [iter delta])");
		return;
	}
	//Points
	vector<CParticleF> points;
	const int* dimsP;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for(int i=0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			points.push_back(CParticleF(x, y));
		}
	}
	int K = 1;
	if(nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(K,prhs[1],classMode);
	} 

	points = removeDegenerate(points);

	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	vector<MovingParticle*> particles = createMovingParticles(points);
	map<MovingParticle*,Vertex<graphKey>*> hmap;
	for(int i=0; i<particles.size(); ++i)
	{
		printf("%d: (%f,%f), (%f,%f), %s\n", 
			i+1, particles[i]->p.m_X, particles[i]->p.m_Y, particles[i]->v.m_X, particles[i]->v.m_Y, particles[i]->convex ? "convex": "reflex");
		Vertex<graphKey>* vt = factory->makeVertex(points[i]);
		hmap[particles[i]] = vt;
	}
	for(int j=0; j<particles.size(); ++j)
	{
		int j2 = (j+1) % particles.size();
		factory->makeEdge(hmap[particles[j]], hmap[particles[j2]], 0.0f);
	}

	StraightAxisPropagation(particles, hmap);

	vector<Vertex<graphKey>*> vertices = factory->vertices;
	vector<Edge<graphKey>*> edges = factory->edges;
	map<Vertex<graphKey>*,int> vmap;
	for(int i=0; i<vertices.size(); ++i)
	{
		vmap[vertices[i]] = i;
	}

	if(nlhs >= 1)
	{
		const int dims[] = {vertices.size(), 2};
		vector<float> F(dims[0]*dims[1]);
		for(int i=0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], vertices[i]->key.m_X);
			SetData2(F, i, 1, dims[0], dims[1], vertices[i]->key.m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if(nlhs >= 2)
	{
		const int dims[] = {edges.size(), 3};
		vector<int> F(dims[0]*dims[1], 0);
		for(int i=0; i<edges.size(); ++i)
		{
			int j1 = vmap[edges[i]->u];
			int j2 = vmap[edges[i]->v];
			SetData2(F, i, 0, dims[0], dims[1], (j1+1));
			SetData2(F, i, 1, dims[0], dims[1], (j2+1));
			SetData2(F, i, 2, dims[0], dims[1], (int)(edges[i]->w));
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	factory->Clean();
	mexUnlock();
}

