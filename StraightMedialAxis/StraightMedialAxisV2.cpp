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
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <TriangulationWithMemento.h>
#include <FragmentInfo.h>
//#include <TriangulationHelper.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
typedef CParticleF graphKey;

GraphFactory<graphKey>* GraphFactory<graphKey>::_instance = NULL;

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

struct MovingParticle
{
	MovingParticle(){}
	MovingParticle(vector<CParticleF>& points, int k)
	{
		this->p = points[k];
		int k0 = (k-1+points.size())  % points.size();
		int k2 = (k+1) % points.size();
		
		CParticleF m((points[k0].m_X+points[k2].m_X)/2, (points[k0].m_Y+points[k2].m_Y)/2);
		this->convex = inside(m, points);
		this->v = _calculateVelocity(p, points[k0], points[k2]);
	}
	CParticleF move(float t)
	{
		p.m_X = p.m_X + t * v.m_X;
		p.m_Y = p.m_Y + t * v.m_Y;
		return p;
	}
	CParticleF _calculateVelocity(CParticleF& o, CParticleF& a, CParticleF& b)
	{
		CParticleF x = NormalizedDirection(a, o);
		CParticleF y = NormalizedDirection(b, o);
		CParticleF z((x.m_X+y.m_X)/2, (x.m_Y+y.m_Y)/2);
		float vx = z.m_X;
		float vy = z.m_Y;
		float len0 = sqrt(vx*vx + vy*vy);
		vx = vx / len0;
		vy = vy / len0;

		float ang = (PI - GetVisualAngle(a.m_X, a.m_Y, b.m_X, b.m_Y, o.m_X, o.m_Y)) / 2.0;
		float len = 1.0f / cos(ang);
		vx = vx * len;
		vy = vy * len;
		if(convex == false)
		{
			vx = -vx;
			vy = -vy;
		}
		return CParticleF(vx, vy);
	}

	CParticleF p;
	CParticleF v;
	bool convex;
};

struct EventDescription
{
	float t;
	int p; //an index to the even particle.
	int q; //an index to a colliding particle if edge event. the first particle of the colliding segment if split event.
	CParticleF loc; //the location of the event
	//bool edgeEvent;
};

bool
coLinear(CParticleF& a, CParticleF& b, CParticleF& c, float precision = 1.0e-3)
{
	float d = Distance2Line(a, c, b);
	return d < precision;
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

EventDescription
nextEdgeEvent(vector<MovingParticle>& v)
{
	EventDescription ev;
	ev.t = std::numeric_limits<float>::infinity();
	for(int i=0; i<v.size(); ++i)
	{
		if(v[i].convex == false) continue;
		CParticleF qi(v[i].p.m_X+v[i].v.m_X, v[i].p.m_Y+v[i].v.m_Y);
		for(int j=i+1; j<v.size(); ++j)
		{
			if(v[j].convex == false) continue;
			CParticleF qj(v[j].p.m_X+v[j].v.m_X, v[j].p.m_Y+v[j].v.m_Y);
			pair<float,float> param = _IntersectConvexPolygon::intersect(v[i].p, qi, v[j].p, qj);
			float t = Max(param.first, param.second);
			if(t>=0 && t < ev.t)
			{
				ev.t = t;
				ev.p = i;
				ev.q = j;
				ev.loc = CParticleF((1-param.first)*v[i].p.m_X+param.first*qi.m_X, (1-param.first)*v[i].p.m_Y+param.first*qi.m_Y); 
			}
		}
	}
	return ev;
}

/*
For the event to occur, a non-convex particle touches moving segment that is not adjacent to the particle.
To detect such event, we consider the moving line segment as a plane with the third axis as time.
We then compute the time the particle (in 3D) intersects the plane.
We then need to check if the time is valid. It has to be positive and has to come after a particle intersects
with traces of both particles of the line segment.
*/
EventDescription
nextSplitEvent(vector<MovingParticle>& v)
{
	EventDescription ev;
	ev.t = std::numeric_limits<float>::infinity();
	for(int i=0; i<v.size(); ++i)
	{
		if(v[i].convex) continue;
		CParticleF ol = v[i].p;
		CParticleF ul = v[i].v;
		ul.m_Z = 1.0;
		for(int j=0; j<v.size(); ++j)
		{
			int j0 = (j-1+v.size()) % v.size();
			int j2 = (j+1) % v.size();
			if(i==j || i==j2 || i==j0) continue;
			CParticleF op = v[j].p;
			CParticleF oq = v[j2].p;
			CParticleF ov = NormalizedDirection(op, oq);
			CParticleF up0(ov.m_Y, -ov.m_X, 1.0f);
			CParticleF up = crossProduct(ov, up0);
			if(up.m_X*v[j].v.m_X + up.m_Y*v[j].v.m_Y < 0)
			{
				up.m_X = -up.m_X;
				up.m_Y = -up.m_Y;
			}
			float t = intersectPlaneAndLine(op, up, ol, ul);
			if(i+1==10 && j+1==1)
			{
				printf("Plane: (%f,%f,%f) (%f,%f,%f).\n", op.m_X, op.m_Y, op.m_Z, up.m_X, up.m_Y, up.m_Z);
				printf("Line: (%f,%f,%f), (%f,%f,%f).\n", ol.m_X, ol.m_Y, ol.m_Z, ul.m_X, ul.m_Y, ul.m_Z);
				printf("t=%f, location=(%f,%f,%f).\n", t, ol.m_X+t*ul.m_X, ol.m_Y+t*ul.m_Y, ol.m_Z+t*ul.m_Z);
			}
			pair<float,float> st0 = _IntersectConvexPolygon::intersect(ol, CParticleF(ol.m_X+ul.m_X, ol.m_Y+ul.m_Y), 
				op, CParticleF(op.m_X+v[j].v.m_X, op.m_Y+v[j].v.m_Y));
			pair<float,float> st2 = _IntersectConvexPolygon::intersect(ol, CParticleF(ol.m_X+ul.m_X, ol.m_Y+ul.m_Y), 
				oq, CParticleF(oq.m_X+v[j2].v.m_X, oq.m_Y+v[j2].v.m_Y));
			float t0=st0.first;
			float t2=st2.first;
			if(t >= 0 && (t>t0 && t>t2) && t < ev.t)
			{
				ev.t = t;
				ev.p = i;
				ev.q = j;
				ev.loc = CParticleF(ol.m_X + t*ul.m_X, ol.m_Y + t*ul.m_Y); 
			}
		}
	}
	return ev;
}

vector<MovingParticle>
createMovingParticles(vector<CParticleF>& points)
{
	vector<MovingParticle> updated;
	for(int i=0; i<points.size(); ++i)
	{
		updated.push_back(MovingParticle(points, i));
	}
	return updated;
}

vector<MovingParticle>
updateEdgeEvent(vector<MovingParticle>& v, 
				EventDescription& ev)
{
	int k = ev.p;
	vector<CParticleF> points;
	for(int i=0; i<v.size(); ++i)
	{
		if(i==k) continue;
		points.push_back(v[i].p);
	}
	return createMovingParticles(points);
}

vector<int> 
getFirstSplitIndices(EventDescription& ev, int size)
{
	vector<int> ind;
	int i=ev.p;
	int end = ev.q;
	do
	{
		ind.push_back(i);
		i = (i+1) % size;
	}
	while(i != end);
	return ind;
}

vector<int> 
getSecondSplitIndices(EventDescription& ev, int size)
{
	vector<int> ind;
	int i=ev.p;
	int end = (ev.q+1) % size;
	do
	{
		ind.push_back(i);
		i = (i+1) % size;
	}
	while(i != end);
	return ind;
}

pair<vector<MovingParticle>,vector<MovingParticle>>
updateSplitEvent(vector<MovingParticle>& v, 
				EventDescription& ev)
{
	vector<int> ind1 = getFirstSplitIndices(ev, v.size());
	vector<CParticleF> points1;
	for(int i=0; i<ind1.size(); ++i)
	{
		points1.push_back(v[ind1[i]].p);
	}

	vector<int> ind2 = getSecondSplitIndices(ev, v.size());
	vector<CParticleF> points2;
	for(int i=0; i<ind2.size(); ++i)
	{
		points2.push_back(v[ind2[i]].p);
	}

	pair<vector<MovingParticle>,vector<MovingParticle>> splitted;
	splitted.first = createMovingParticles(points1);
	splitted.second = createMovingParticles(points2);
	return splitted;
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

	vector<vector<MovingParticle>> shapes;
	vector<vector<Vertex<graphKey>*>> vshapes;
	GraphFactory<graphKey>* factory = GraphFactory<graphKey>::GetInstance();
	vector<Vertex<graphKey>*> vfront;
	vector<MovingParticle> particles;
	for(int i=0; i<points.size(); ++i)
	{
		particles.push_back(MovingParticle(points, i));
		Vertex<graphKey>* v = factory->makeVertex(points[i]);
		vfront.push_back(v);
	}
	for(int i=0; i<vfront.size(); ++i)
	{
		vfront[i]->Add(factory->makeEdge(vfront[i], vfront[(i+1)%vfront.size()], 0.0f));
	}
	shapes.push_back(particles);
	vshapes.push_back(vfront);

	for(int k=1; k<=K; ++k)
	{
		int sz = shapes.size();
		vector<pair<CParticleF,CParticleF>> trace;
		for(int j=0; j<sz; ++j)
		{
			EventDescription ev1 = nextEdgeEvent(shapes[j]);
			EventDescription ev2 = nextSplitEvent(shapes[j]);
			float t = Min(ev1.t, ev2.t);
			if(t < std::numeric_limits<float>::infinity())
			{
				for(int i=0; i<shapes[j].size(); ++i)
				{
					pair<CParticleF,CParticleF> p;
					p.first = shapes[j][i].p;
					shapes[j][i].move(t);
					p.second = shapes[j][i].p;
					trace.push_back(p);
				}
				if(ev1.t < ev2.t)
				{
					printf("merging points %d and %d at time %f at location (%f,%f).\n", ev1.p+1, ev1.q+1, ev1.t, ev1.loc.m_X, ev1.loc.m_Y);
					shapes[j] = updateEdgeEvent(shapes[j], ev1);
				}
				else
				{
					printf("Splitting by point %d, segment %d-%d at time %f at location (%f,%f).\n", 
						ev2.p+1, ev2.q+1, ((ev2.q+1) % shapes[j].size()) + 1, ev2.t, ev2.loc.m_X, ev2.loc.m_Y);
					pair<vector<MovingParticle>,vector<MovingParticle>> p = updateSplitEvent(shapes[j], ev2);
					shapes[j] = p.first;
					shapes.push_back(p.second);
				}
			}
		}
		vector<Vertex<graphKey>*> vfront;
		vector<CParticleF> dest;
		for(int j=0; j<trace.size(); ++j)
		{
			CParticleF p = trace[j].first;
			CParticleF q = trace[j].second;
			int k1 = distance(points.begin(), find(points.begin(), points.end(), p));
			int k2 = distance(dest.begin(), find(dest.begin(), dest.end(), q));
			if(k2==dest.size())
			{
				Vertex<graphKey>* vertex = factory->makeVertex(q);
				vfront.push_back(vertex);
				dest.push_back(q);
			}
			vshapes[k-1][k1]->Add(factory->makeEdge(vshapes[k-1][k1], vfront[k2]));
		}
		vshapes.push_back(vfront);
		points = dest;
	}

	vector<Vertex<graphKey>*> vertices;
	vector<Edge<graphKey>*> edges;
	map<Vertex<graphKey>*,int> vmap;
	for(int i=0; i<vshapes.size(); ++i)
	{
		for(int j=0; j<vshapes[i].size(); ++j)
		{
			vertices.push_back(vshapes[i][j]);
			vmap[vshapes[i][j]] = vertices.size()-1;
			for(int k=0; k<vshapes[i][j]->aList.size(); ++k)
			{
				edges.push_back(vshapes[i][j]->aList[k]);
			}
		}
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
		const int dims[] = {edges.size(), 2};
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

