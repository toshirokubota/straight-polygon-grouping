#include <set>
#include <map>
#include <TriangulationSaliencyUtility.h>
#include <szMiscOperations.h>
#include <TriangulationHelper.h>
#include <szMexUtility.h>
#include <Eigen.h>

float lengthTerm(float length, float lengthThres, float alpha, float angleTerm)
{
	//lengthThres *= 1.0 + exp(-alpha*(1.0-angleTerm));
	return 1.0-1/(1+exp(-alpha*(length-lengthThres)/lengthThres));
}

float angleTerm(float angle)
{
	//return (1 - cos(angle))*(1 - cos(angle))/4.0f;
	float at = (1 - cos(angle))/2.0f;
	return at; // at * at; // * at * at * at;
}

float angleTerm(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2)
{
	return angleTerm(angleBetween(e1, e2));
}

float 
pairFitness(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2, float beta, float gamma)
{
	float a1 = angleTerm(e1, e2);
	float len1 = lengthTerm(e1->Length(), beta, gamma, a1);
	float len2 = lengthTerm(e2->Length(), beta, gamma, a1);
	return len1 * len2 * a1;
}

vector<Triangulation::_Internal::_vertex*>
orientVertices(vector<Triangulation::_Internal::_vertex*>& vertices, CParticleF& o)
{
	vector<pair<float,Triangulation::_Internal::_vertex*>> pairs;
	for(int i=0; i<vertices.size(); ++i)
	{
		float angle = GetVisualDirection(vertices[i]->p.m_X, vertices[i]->p.m_Y, o.m_X, o.m_Y);
		pairs.push_back(pair<float,Triangulation::_Internal::_vertex*>(angle, vertices[i]));
	}
	sort(pairs.begin(), pairs.end());
	vector<Triangulation::_Internal::_vertex*> arranged;
	for(int i=0; i<vertices.size(); ++i)
	{
		if(i>0)
		{
			if(vertices[i] == vertices[i-1]) continue; //remove duplicates
		}
		arranged.push_back(pairs[i].second);
	}
	return arranged;
}

bool
isBoundary(Triangulation::_Internal::_vertex* u)
{
	for(int i=0; i<u->edges.size(); ++i)
	{
		if(u->edges[i]->type == Triangulation::_Internal::Boundary)
		{
			return true;
		}
	}
	return false;
}

float
pointFitness(Triangulation::_Internal::_vertex* u, 
			 Triangulation::_Internal::_vertex* v,
			 Triangulation::_Internal::_vertex* w,
			 float beta, float gamma)
{
	if((int)u->p.m_X==91 && (int)u->p.m_Y==6)
	{
		u->p.m_Life += 0;
	}
	Triangulation::_Internal::_edge* e1 = findEdge(u, v);
	Triangulation::_Internal::_edge* e2 = findEdge(u, w);
	if(e1 == NULL || e2==NULL) return 0.0f;
	
	set<Triangulation::_Internal::_vertex*> vset;
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = u->edges[i];
		vset.insert(e->vertices[0] == u ? e->vertices[1]: e->vertices[0]);
	}
	for(int i=0; i<v->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = v->edges[i];
		vset.insert(e->vertices[0] == v ? e->vertices[1]: e->vertices[0]);
	}
	for(int i=0; i<w->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = w->edges[i];
		vset.insert(e->vertices[0] == w ? e->vertices[1]: e->vertices[0]);
	}
	vector<Triangulation::_Internal::_vertex*> vertices;
	for(set<Triangulation::_Internal::_vertex*>::iterator p=vset.begin(); p!=vset.end(); ++p)
	{
		vertices.push_back(*p);
	}
	vertices = orientVertices(vertices, u->p);
	vector<CParticleF> left;
	{
		int k = distance(vertices.begin(), find(vertices.begin(), vertices.end(), v));
		k = (k + 1) % vertices.size();
		while(vertices[k] != w)
		{
			left.push_back(vertices[k]->p);
			k = (k + 1) % vertices.size();
		}
	}
	vector<CParticleF> right;
	{
		int k = distance(vertices.begin(), find(vertices.begin(), vertices.end(), w));
		k = (k + 1) % vertices.size();
		while(vertices[k] != v)
		{
			right.push_back(vertices[k]->p);
			k = (k + 1) % vertices.size();
		}
	}
	float leftTerm;
	float rightTerm;
	if(left.empty())
	{
		//if this is because the triangle is right at the boundary, then return 1.0f.
		if(isBoundary(u))
		{
			leftTerm = 1.0;
		}
		//else there must be an edge between u and w. Find the distance from u to the edge and use the twice the length for the missing side.
		else
		{
			leftTerm = lengthTerm(2.0 * Distance2Line(v->p, w->p, u->p), beta, gamma);
		}
	}
	else
	{
		float sum = 0;
		for(int i=0; i<left.size(); ++i)
		{
			sum += Distance(u->p, left[i]);
		}
		float av = sum / left.size();
		leftTerm = 1.0 - lengthTerm(av, beta, gamma);
	}
	if(right.empty())
	{
		//if this is because the triangle is right at the boundary, then return 1.0f.
		if(isBoundary(u))
		{
			rightTerm = 1.0f;
		}
		else
		{
			rightTerm = lengthTerm(2.0 * Distance2Line(v->p, w->p, u->p), beta, gamma);
		}
	}
	else
	{
		float sum = 0;
		for(int i=0; i<right.size(); ++i)
		{
			sum += Distance(u->p, right[i]);
		}
		float av = sum / right.size();
		rightTerm = 1.0 - lengthTerm(av, beta, gamma);
	}
	return Max(leftTerm, rightTerm);
	//return exp(-4.0 * (1.0-leftTerm) * (1.0-rightTerm));
}

/*
Perform 2D eigenvalue/eigenvector computation from a set of points with a weight.
The weight of the point is stored in m_Life.
The first two components are eigenvalues (largest then smallest).
The next two components are eigenvectors of the largest eigenvalue.
The next two components are eigenvectors of the smallest eigenvalue.
The next two components are the center of the point cloud.
*/
vector<float>
Moment(vector<CParticleF>& pnts)
{
	float x=0, y=0, xx=0, yy=0, xy=0;
	for(int i=0; i<pnts.size(); ++i)
	{
		x += pnts[i].m_X * pnts[i].m_Life;
		y += pnts[i].m_Y * pnts[i].m_Life;
		xx += pnts[i].m_X * pnts[i].m_X * pnts[i].m_Life;
		yy += pnts[i].m_Y * pnts[i].m_Y * pnts[i].m_Life;
		xy += pnts[i].m_X * pnts[i].m_Y * pnts[i].m_Life;
	}
	x /= pnts.size();
	y /= pnts.size();
	xx = xx / pnts.size() - x * x;
	yy = yy / pnts.size() - y * y;
	xy = xy / pnts.size() - x * y;

	float det = sqrt((xx-yy)*(xx-yy) + 4*xy*xy);
	float e1 = (xx + yy + det)/2.0;
	float e2 = (xx + yy - det)/2.0;
	vector<float> m;
	m.push_back(e1);
	m.push_back(e2);
	{

		float xv = -(-xx + yy - det)/(2.*xy);
		float len = sqrt(xv*xv + 1);
		m.push_back(xv/len);
		m.push_back(1.0/len);
		float xv2 = -(-xx + yy + det)/(2.*xy);
		float len2 = sqrt(xv2*xv2 + 1);
		m.push_back(xv2/len2);
		m.push_back(1.0/len2);
	}
	m.push_back(x);
	m.push_back(y);
	return m;
}

/*
Similar to Moment(...) defined above.
This one does not use m_Life, only m_X and m_Y.
*/
vector<float>
Moment2(vector<CParticleF>& pnts)
{
	float x=0, y=0, xx=0, yy=0, xy=0;
	for(int i=0; i<pnts.size(); ++i)
	{
		x += pnts[i].m_X;
		y += pnts[i].m_Y;
	}
	x /= pnts.size();
	y /= pnts.size();
	for(int i=0; i<pnts.size(); ++i)
	{
		xx += (pnts[i].m_X-x) * (pnts[i].m_X-x);
		yy += (pnts[i].m_Y-y) * (pnts[i].m_Y-y);
		xy += (pnts[i].m_X-x) * (pnts[i].m_Y-y);
	}
	xx = xx / pnts.size();// - x * x;
	yy = yy / pnts.size();// - y * y;
	xy = xy / pnts.size();// - x * y;*/
	//printf("%f %f %f %f %f\n", x, y, xx, yy, xy);

	float det = sqrt((xx-yy)*(xx-yy) + 4*xy*xy);
	float e1 = (xx + yy + det)/2.0;
	float e2 = (xx + yy - det)/2.0;
	vector<float> m;
	m.push_back(e1);
	m.push_back(e2);
	if(Abs(xy) < 1.0e-8){
		m.push_back(1);
		m.push_back(0);
		m.push_back(0);
		m.push_back(1);
	}
	else
	{
		float xv = -(-xx + yy - det)/(2*xy); 
		float len = sqrt(xv*xv + 1);
		m.push_back(xv/len);
		m.push_back(1.0/len);
		float xv2 = -(-xx + yy + det)/(2*xy);
		float len2 = sqrt(xv2*xv2 + 1);
		m.push_back(xv2/len2);
		m.push_back(1.0/len2);
	}
	m.push_back(x);
	m.push_back(y);
	return m;
}

/*
Circular fit of points.
[0] = radius.
[1], [2]: center
*/
vector<float>
MomentCircular(vector<CParticleF>& pnts)
{
	float x=0, y=0, xx=0, yy=0, xy=0;
	for(int i=0; i<pnts.size(); ++i)
	{	
		x += pnts[i].m_X;
		y += pnts[i].m_Y;
		xx += pnts[i].m_X * pnts[i].m_X;
		yy += pnts[i].m_Y * pnts[i].m_Y;
	}
	x /= pnts.size();
	y /= pnts.size();
	xx = xx / pnts.size() - x * x;
	yy = yy / pnts.size() - y * y;

	vector<float> m;
	m.push_back(sqrt(xx + yy));
	m.push_back(x);
	m.push_back(y);
	return m;
}

float
pointRadialFitness(Triangulation::_Internal::_vertex* u, 
					float beta, float gamma)
{
	float eps = 1.0e-6;
	map<Triangulation::_Internal::_vertex*,float> vmap;
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = u->edges[i];
		Triangulation::_Internal::_vertex* v = e->vertices[0] == u ? e->vertices[1]: e->vertices[0];
		vmap[v] = 0;
	}
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e1 = u->edges[i];
		Triangulation::_Internal::_vertex* v = e1->vertices[0] == u ? e1->vertices[1]: e1->vertices[0];
		for(int j=i+1; j<u->edges.size(); ++j)
		{
			Triangulation::_Internal::_edge* e2 = u->edges[j];
			Triangulation::_Internal::_vertex* w = e2->vertices[0] == u ? e2->vertices[1]: e2->vertices[0];
			float fit = pairFitness(e1, e2, beta, gamma);
			vmap[v] = vmap[v] + fit;
			vmap[w] = vmap[w] + fit;
		}
	}
	vector<CParticleF> pnts;
	for(int i=0; i<u->edges.size(); ++i)
	{
		Triangulation::_Internal::_edge* e = u->edges[i];
		Triangulation::_Internal::_vertex* v = e->vertices[0] == u ? e->vertices[1]: e->vertices[0];
		float theta = atan2(v->p.m_Y-u->p.m_Y, v->p.m_X-u->p.m_X);
		pnts.push_back(CParticleF(cos(theta), sin(theta), 0, vmap[v]));
		//pnts.push_back(CParticleF(v->p.m_X-u->p.m_X, v->p.m_Y-u->p.m_Y, 0, vmap[v]));
	}
	vector<float> m = Moment(pnts);
	return 1.0 - m[1] / Max(eps, m[0]);
}

/*StructuralTensor3D::StructuralTensor3D(vector<CParticleF>& pnts)
{
	if(pnts.empty()) 
	{
		_init();
		return;
	}

	double x=0, y=0, z=0, xx=0, yy=0, zz=0, xy=0, yz=0, xz=0;
	for(int i=0; i<pnts.size(); ++i)
	{
		float x0=pnts[i].m_X;
		float y0=pnts[i].m_Y;
		float z0=pnts[i].m_Z;
		x+=x0;
		y+=y0;
		z+=z0;
		xx+=x0*x0;
		yy+=y0*y0;
		zz+=z0*z0;
		xy+=x0*y0;
		yz+=y0*z0;
		xz+=x0*z0;
	}
	x /= pnts.size();
	y /= pnts.size();
	z /= pnts.size();
	xx = xx/pnts.size() - x*x;
	yy = yy/pnts.size() - y*y;
	zz = zz/pnts.size() - z*z;
	xy = xy/pnts.size() - x*y;
	yz = yz/pnts.size() - y*z;
	xz = xz/pnts.size() - x*z;
	vector<double> M(9);
	M[0]=xx; M[1]=xy; M[2]=xz;
	M[3]=xy; M[4]=yy; M[5]=yz;
	M[6]=xz; M[7]=yz; M[8]=zz;

	struct_eigen st = computeEigenValues(M, 3);
	this->cx = x;
	this->cy = y;
	this->cz = z;
	for(int i=0; i<3; ++i)
	{
		this->evals[i] = st.evals[i];
		for(int j=0; j<3; ++j)
		{
			this->axes[i][j] = st.evecs[i][j];
		}
	}
}*/
