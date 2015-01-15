#include <intersectionConvexPolygons.h>
#include <szMiscOperations.h>
//#include <DotGroupUtility.h>
#include <szMexUtility.h>
#include <szConvexHull2D.h>
#include <cassert>
namespace _IntersectConvexPolygon
{
	/*
	*/
	/*bool
	inside(const CParticleF& p, const vector<CParticleF>& pnts)
	{
		bool c = false;
		int x = p.m_X;
		int y = p.m_Y;
		for (int i = 0, j = pnts.size()-1; i < pnts.size(); j = i++) 
		{
			float xi = pnts[i].m_X;
			float yi = pnts[i].m_Y;
			float xj = pnts[j].m_X;
			float yj = pnts[j].m_Y;
			if ( ((yi>y) != (yj>y)) &&
				(x < (xj-xi) * (y-yi) / (yj-yi) + xi) )
				c = !c;
		}
		return c;
	}*/

	/*
	return true if p is on the boundary of a shape delineated by Q
	*/
	bool
	onBoundary(const CParticleF& p, const std::vector<CParticleF>& Q, float eps)
	{
		for (int i = 0; i < Q.size(); ++i)
		{
			float d = Distance2LineSegment(Q[i], Q[i == Q.size() - 1 ? 0 : i + 1], p);
			if (d < eps)
			{
				return true;
			}
		}
		return false;
	}

	/*
	find an intersection between two line segments defined by [p1 - p2] and [q1 - q2].
	It returns a parameter s and t such that the intersection point is (1-s)*p1 + s*p2
	and (1-t)*q1+t*q2.
	*/
	pair<float,float>
	intersect(const CParticleF& p1, const CParticleF& p2, const CParticleF& q1, const CParticleF& q2)
	{
		double x0=p1.m_X;
		double y0=p1.m_Y;
		double x1=p2.m_X;
		double y1=p2.m_Y;
		double x2=q1.m_X;
		double y2=q1.m_Y;
		double x3=q2.m_X;
		double y3=q2.m_Y;
		double det = (x0-x1)*(y2-y3) - (x2-x3)*(y0-y1);

		double s = (x0*(y2-y3)+x2*(y3-y0)+x3*(y0-y2))/det;
		double t = (x0*(y2-y1)+x1*(y0-y2)+x2*(y1-y0))/det;

		pair<float,float> param((float)s,(float)t);
		return param;
	}

	/*
	get a bonding box as a vector of CParticleF.
	*/
	/*std::vector<CParticleF>
	boundingBox(std::vector<CParticleF>& P)
	{
		std::vector<CParticleF> box(4);
		if (P.empty())
		{
			for (int i = 0; i < 4; ++i)
			{
				box[i].m_X = std::numeric_limits<float>::quiet_NaN;
				box[i].m_Y = std::numeric_limits<float>::quiet_NaN;
			}
		}
		else
		{
			float minxp = P[0].m_X, minyp = P[0].m_Y;
			float maxxp = P[0].m_X, maxyp = P[0].m_Y;
			for (int i = 1; i < P.size(); ++i)
			{
				minxp = Min(minxp, P[i].m_X);
				minyp = Min(minyp, P[i].m_Y);
				maxxp = Max(maxxp, P[i].m_X);
				maxyp = Max(maxyp, P[i].m_Y);
			}
			box[0] = CParticleF(minxp, minyp);
			box[1] = CParticleF(minxp, maxyp);
			box[2] = CParticleF(maxxp, maxyp);
			box[3] = CParticleF(maxxp, minyp);
		}
		return box;
	}*/

	/*
	return true if bounding boxes of P and Q do not overlap.
	*/
	/*bool
	intersectingBounds(std::vector<CParticleF>& P, 
						 std::vector<CParticleF>& Q)
	{
		if (P.empty() || Q.empty())
		{
			return true;
		}
		vector<CParticleF> boxP = boundingBox(P);
		vector<CParticleF> boxQ = boundingBox(Q);
		float minxp = boxP[0].m_X, minyp = boxP[0].m_Y;
		float maxxp = boxP[2].m_X, maxyp = boxP[2].m_Y;
		float minxq = boxQ[0].m_X, minyq = boxQ[0].m_Y;
		float maxxq = boxQ[2].m_X, maxyq = boxQ[2].m_Y;
		if (minxp > maxxq) return false;
		else if (maxxp < minxq) return false;
		else if (minyp > maxyq) return false;
		else if (maxyp < minyq) return false;
		else return true;
	}*/

	bool
	containsTheOther(const std::vector<CParticleF>& P, 
						 const std::vector<CParticleF>& Q)
	{
		for(int i=0; i<Q.size(); ++i)
		{
			if(inside(Q[i], P)==false) return false;
		}
		return true;
	}

	/*
	a brute force approach to finding the intersection of two convex polygons.
	1. add all points in P that are inside Q to a set.
	2. add all points in Q that are inside P to the set.
	3. add all intersection points to the set.
	4. return the convex hull of the set.
	*/
	std::vector<CParticleF> 
	intersectConvexHulls(const std::vector<CParticleF>& P, 
						 const std::vector<CParticleF>& Q)
	{
		/*if (intersectingBounds(P, Q) == false)
		{
			return std::vector<CParticleF>();
		}*/
		if(containsTheOther(P, Q))
		{
			return Q;
		}
		else if(containsTheOther(Q, P))
		{
			return P;
		}
		else
		{
			std::vector<CParticleF> pset;
			for(int i=0; i<P.size(); ++i)
			{
				if(inside(P[i],Q) || onBoundary(P[i], Q, 1.0e-5))
				{
					pset.push_back(P[i]);
				}
			}
			for(int i=0; i<Q.size(); ++i)
			{
				if (inside(Q[i], P) || onBoundary(Q[i], P, 1.0e-5))
				{
					pset.push_back(Q[i]);
				}
			}
			for(int i=0; i<P.size(); ++i)
			{
				int j = (i+1) % P.size();
				if(inside(P[i], Q)==false || inside(P[j],Q)==false)
				{
					for(int i2=0; i2<Q.size(); ++i2)
					{
						int j2 = (i2+1) % Q.size();
						pair<float,float> param = intersect(P[i], P[j], Q[i2], Q[j2]);
						if(param.first>0 && param.first<1.0f && param.second>0 && param.second<1.0f)
						{
							float t = param.first;
							CParticleF p((1-t)*P[i].m_X+t*P[j].m_X, (1-t)*P[i].m_Y+t*P[j].m_Y);
							pset.push_back(p);
						}
					}
				}
			}
			return ConvexHull2D(pset);
		}
	}
}
