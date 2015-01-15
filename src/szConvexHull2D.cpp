#include <szConvexHull2D.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMiscOperations.h>

CParticleF
FindRightMostLowestPoint(const vector<CParticleF>& vp)
{
	CParticleF cand = vp[0];
	for(int i=1; i<vp.size(); ++i)
	{
		if(cand.m_Y > vp[i].m_Y || (cand.m_Y == vp[i].m_Y && cand.m_X < vp[i].m_X))
		{
			cand = vp[i];
		}
	}
	return cand;
}

float
isLeft( const CParticleF& P0, const CParticleF& P1, const CParticleF& P2 )
{
	return (P1.m_X - P0.m_X)*(P2.m_Y - P0.m_Y) - (P2.m_X - P0.m_X)*(P1.m_Y - P0.m_Y);
}

float
getAngle(float x, float y)
{
	return atan2(y, x);
}

vector<CParticleF>
ConvexHull2D(const vector<CParticleF>& vp)
{
	if(vp.size() <= 3) return vp;

	CParticleF p0 = FindRightMostLowestPoint(vp);
	vector<CParticleF> vfront = vp;
	vector<CParticleF>::iterator ip = find(vfront.begin(), vfront.end(), p0);
	vfront.erase(ip);

	//sort particles based on its angle from the right-most-lowest point
	for(int i=0; i<vfront.size(); ++i)
	{
		float angle = getAngle(vfront[i].m_X-p0.m_X, vfront[i].m_Y-p0.m_Y);
		vfront[i].m_Life = angle;
	}
	sort(vfront.begin(), vfront.end());

	//remove tie
	vector<CParticleF> vp2;
	int n=0;
	while(n<vfront.size())
	{
		CParticleF p = vfront[n];
		int m=n+1;
		while(m<vfront.size() && Abs(vfront[m].m_Life - p.m_Life)<=std::numeric_limits<float>::epsilon())
		{
			float d1 = Distance(p, p0);
			float d2 = Distance(vfront[m], p0);
			if(d2 > d1)
			{
				p = vfront[m];
			}
			m++;
		}
		vp2.push_back(p);
		n = m;
	}
	//select those that form convex hull
	vector<CParticleF> hull(vp2.size()+1);
	hull[0] = p0;
	hull[1] = vp2[0];
	int count = 2;
	for(int i=1; i<vp2.size(); ++i)
	{
		while(1)
		{
			CParticleF p1 = hull[count-2];
			CParticleF p2 = hull[count-1];
			if(isLeft(p1, p2, vp2[i])<0)
			{
				count--;
			}
			else
			{
				break;
			}
		}
		hull[count] = vp2[i];
		count++;
	}
	hull.erase(hull.begin()+count, hull.end());
	return hull;
}

void
PaintInsideHull(vector<unsigned char>& im, 
				const vector<CParticleF>& vch, 
				const int* dims)
{
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			bool bOut = false;
			CParticleF q(x, y, 0);
			for(int n=0; n<vch.size(); ++n)
			{
				CParticleF p0 = vch[n];
				CParticleF p1 = vch[(n+1) % vch.size()];
				if(isLeft(p0, p1, q)<0)
				{
					bOut = true;
					break;
				}
			}
			if(bOut == false)
			{
				SetData2(im, x, y, dims[0], dims[1], (unsigned char)1);
			}
		}
	}
}

/*
Check if three 3d points are colinera. I've been a slack and use the difference in the
directional unit vectors as the indicator.
*/
float
_colinear(const CParticleF& p, const CParticleF& q, const CParticleF& r)
{
	CParticleF u = UnitVector(Difference(q, p));
	CParticleF v = UnitVector(Difference(r, p));
	return Distance(u, v);
}

/*
By running this, we can remove points that are 'nearly' co-linear.
The larger the threshold, the more points are considered 'co-linear'.
Empirically speaking, 0.01 seems to work reasonably.
*/
vector<CParticleF>
RemoveColinearPoints(const vector<CParticleF>& hull, float thres)
{
	vector<CParticleF> points = hull;
	while(true)
	{
		bool bChanged = false;
		for(int i=1; i<points.size()-1; ++i)
		{
			if(_colinear(points[i-1], points[i], points[i+1]) < thres)
			{
				points.erase(points.begin()+i);
				bChanged = true;
				break;
			}
		}

		if(bChanged == false) break;
	}
	return points;
}
