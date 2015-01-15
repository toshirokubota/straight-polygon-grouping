#include <MiscGeometry.h>
#include <szMiscOperations.h>
#include <szmexutilitytemplate.h>
#include <MovingParticle.h>
/*
A unit vector from o to x in 2D
*/
CParticleF
NormalizedDirection(CParticleF& x, CParticleF& o)
{
	float dx = x.m_X - o.m_X;
	float dy = x.m_Y - o.m_Y;
	float len = sqrt(dx*dx + dy*dy);
	if (len > 1.0e-5)
	{
		return CParticleF(dx / len, dy / len);
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
	if (len > 1.0e-5)
	{
		return CParticleF(dx / len, dy / len, dz / len);
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
	CParticleF z((x.m_X + y.m_X) / 2, (x.m_Y + y.m_Y) / 2);
	float vx = z.m_X;
	float vy = z.m_Y;
	float len0 = sqrt(vx*vx + vy*vy);
	if (len0 <= 1.0e-5) //this is a colinear point. 
	{
		float ang = GetVisualDirection(b.m_X, b.m_Y, a.m_X, a.m_Y) - PI / 2.0;
		vx = cos(ang);
		vy = sin(ang);
	}
	else
	{
		vx = vx / len0;
		vy = vy / len0;
	}
	CParticleF bs(o.m_X + vx, o.m_Y + vy);

	//There are two bisector directions. 
	//We consistently choose one that is in clock-wise direction.
	vector<CParticleF> pnts(4);
	pnts[0] = o;
	pnts[1] = a;
	pnts[2] = bs;
	pnts[3] = b;
	if (ClockWise(pnts)<0)
	{
		vx = -vx;
		vy = -vy;
	}
	return CParticleF(vx, vy);
}

/*
OP:  a point on the plane.
Z: a perpendicular vector to the plane.
OL: a point on the line.
U: a vector parallel to the line.

It finds a parameter t such that ol + t*u intersects the plane.
*/
float
intersectPlaneAndLine(CParticleF op, CParticleF& z, CParticleF& ol, CParticleF u)
{
	//move the origin to ol.
	op.m_X -= ol.m_X;
	op.m_Y -= ol.m_Y;
	op.m_Z -= ol.m_Z;
	float denom = z.m_X * u.m_X + z.m_Y * u.m_Y + z.m_Z * u.m_Z;
	float numer = z.m_X * op.m_X + z.m_Y * op.m_Y + z.m_Z * op.m_Z;
	//float t = (z.m_X * op.m_X + z.m_Y * op.m_Y + z.m_Z * op.m_Z) / (z.m_X * u.m_X + z.m_Y * u.m_Y + z.m_Z * u.m_Z);
	float t = numer / denom;
	return t;
}

CParticleF
crossProduct(CParticleF& a, CParticleF& b)
{
	return CParticleF(a.m_Y*b.m_Z - a.m_Z*b.m_Y, a.m_Z*b.m_X - a.m_X*b.m_Z, a.m_X*b.m_Y - a.m_Y*b.m_X);
}


float
dotProduct(CParticleF& a, CParticleF& b)
{
	return a.m_X*b.m_X + a.m_Y*b.m_Y + a.m_Z*b.m_Z;
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
	if (z.m_Z < 0)
	{
		z.m_X = -z.m_X;
		z.m_Y = -z.m_Y;
		z.m_Z = -z.m_Z;
	}
	return z;
}

CParticleF
centroidPoint(vector<CParticleF>& poly)
{
	CParticleF c(0, 0, 0);
	for (int i = 0; i<poly.size(); ++i)
	{
		c.m_X += poly[i].m_X;
		c.m_Y += poly[i].m_Y;
		c.m_Z += poly[i].m_Z;
	}
	c.m_X /= (float)poly.size();
	c.m_Y /= (float)poly.size();
	c.m_Z /= (float)poly.size();
	return c;
}

/*
Check if a-b-c lie on a line.
*/
bool
coLinear(CParticleF& a, CParticleF& b, CParticleF& c, float precision)
{
	float d = Distance2Line(a, c, b);
	return d < precision;
}

/*
Check if lines a-b and c-d are parallel.
*/
bool
isParallel(CParticleF& a, CParticleF& b, CParticleF& c, CParticleF& d, float precision)
{
	CParticleF x = UnitVector(CParticleF(b.m_X - a.m_X, b.m_Y - a.m_Y, b.m_Z - a.m_Z));
	CParticleF y = UnitVector(CParticleF(d.m_X - c.m_X, d.m_Y - c.m_Y, d.m_Z - c.m_Z));

	CParticleF z = crossProduct(x, y);
	if (sqrt(z.m_X*z.m_X + z.m_Y*z.m_Y + z.m_Z*z.m_Z) < precision)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*
Check if a moving side (p-q) and a moving particle (a with velocity (vx,vy)) are approaching to each other.
An important assumption is that the side moves to the left when you are at p facing q. This corresponds to q = p->next->p.
*/
bool
isApproaching(CParticleF& a, float vx, float vy, CParticleF& p, CParticleF& q)
{
	//direction of teh moving side
	CParticleF u = perpDirection(p, q);
	CParticleF b = Closest2Line(p, q, a);
	float ip = (a.m_X - b.m_X) * (vx - u.m_X) + (a.m_Y - b.m_Y) * (vy - u.m_Y);
	return ip <= 0; //need the equality not to disregard events that are happening right now
}

/*
Compute a unit vector normal to the side (p->q). The direction is to the left when we are at p and looking toward q.
*/
CParticleF
perpDirection(CParticleF& p, CParticleF& q)
{
	CParticleF nd = NormalizedDirection(q, p);
	return CParticleF(-nd.m_Y, nd.m_X);
}
