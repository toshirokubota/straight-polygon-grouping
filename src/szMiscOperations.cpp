
#include <mex.h>
#include <math.h>
#include <algorithm>
using namespace std;

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>

int GetOffsetX6(int n) 
{
  switch(n) {
  case 1:
    return -1;
  case 2:
    return 1;
  case 0:
  case 3:
  case 4:
  case 5:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetY6(int n) 
{
  switch(n) {
  case 0:
    return -1;
  case 3:
    return 1;
  case 1:
  case 2:
  case 4:
  case 5:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetZ6(int n) 
{
  switch(n) {
  case 4:
    return -1;
  case 5:
    return 1;
  case 0:
  case 1:
  case 2:
  case 3:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
  //const int zoff[6]={0,0,0,0,-1,1};
  //return zoff[n];
}

int GetOffsetX10(int n) 
{
  switch(n) {
  case 1:
  case 4:
  case 6:
    return -1;
  case 3:
  case 5:
  case 8:
    return 1;
  case 0:
  case 2:
  case 7:
  case 9:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetY10(int n) 
{
  switch(n) {
  case 1:
  case 2:
  case 3:
    return -1;
  case 6:
  case 7:
  case 8:
    return 1;
  case 0:
  case 4:
  case 5:
  case 9:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetZ10(int n) 
{
  switch(n) {
  case 0:
    return -1;
  case 9:
    return 1;
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
  //const int zoff[6]={0,0,0,0,-1,1};
  //return zoff[n];
}

int GetOffsetX26(int n) 
{
  if(n<9)
    return (n%3)-1;
  else if(n<17) {
    switch(n) {
    case 9:
    case 12:
    case 14:
      return -1;
    case 10:
    case 15:
      return 0;
    case 11:
    case 13:
    case 16:
      return 1;
	default:
		return 0;
    }
  }
  else 
    return (n-17)%3 - 1;
}
 
int GetOffsetY26(int n) 
{
  if(n<9)
    return (n/3)-1;
  else if(n<17) {
    switch(n) {
    case 9:
    case 10:
    case 11:
      return -1;
    case 12:
    case 13:
      return 0;
    case 14:
    case 15:
    case 16:
      return 1;
	default:
		return 0;
    }
  }
  else 
    return (n-17)/3-1;
}
 
int GetOffsetZ26(int n) {
  if(n<9)
    return -1;
  else if(n<17)
    return 0;
  else //if(n<26)
    return 1;
}

/*
Get the median of the three numbers.
*/
int 
Middle(int a, int b, int c)
{
	if(a >= b && a <=c || a >= c && a <= b) return a;
	else if(b >= a && b <= c || b >= c && b <= a) return b;
	else return c;
}

/*float
Distance(const CParticleF& p1, const CParticleF& p2)
{
	return sqrt((float)(p1.m_X-p2.m_X)*(p1.m_X-p2.m_X)+(p1.m_Y-p2.m_Y)*(p1.m_Y-p2.m_Y));
}*/

float
Distance(float x1, float y1, float x2, float y2)
{
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

/*
Check if a trace of three points (x0,y0) -> (x1,y1) -> (x2,y2) is in Clockwise.
It returns 1 if clockwise, -1 if counter clockwise, and 0 if undetermined.
*/
int
ClockWise(float x0, float y0, float x1, float y1, float x2, float y2)
{
	x1-=x0;
	y1-=y0;
	x2-=x0;
	y2-=y0;
	float cp = x1*y2-x2*y1;
	float eps = 1.0e-6;
	if(cp < -eps)
	{
		return 1;
	}
	else if(cp > eps)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

/*
Check if a trace of points is in Clockwise.
It returns 1 if clockwise, -1 if counter clockwise, and 0 if undetermined.
*/
int
ClockWise(const vector<CParticleF>& points)
{
	float cp = 0;
	for(int i=0, j=points.size()-1; i<points.size(); ++i)
	{
		cp += (points[i].m_X-points[j].m_X)*(points[i].m_Y+points[j].m_Y);
		j=i;
	}

	float eps = 1.0e-6;
	if(cp < -eps)
	{
		return -1;
	}
	else if(cp > eps)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

/*
get the direction of (x, y) as seen from (x0, y0).
*/
float
GetVisualDirection(float x, float y, 
				   float x0, float y0)
{
	float angle = atan2(y-y0, x-x0);
	if(angle != angle) //the contour point is identical to the click point
	{
		angle = 0;
	}
	return angle;
}

/*
Get the angle of view formed by (x1, y1) and (x2, y2) as seen from (x0, y0).
*/
float
GetVisualAngle(float x1, float y1, float x2, float y2, 
			   float x0, float y0)
{
	float a = Distance(x1, y1, x0, y0);
	float b = Distance(x2, y2, x0, y0);
	float c = Distance(x1, y1, x2, y2);
	float eps = 1.0e-8;
	if(a<eps || b<eps)
	{
		return 0;
	}
	else
	{
		float s = (a*a + b*b - c*c)/(2*a*b);
		if (s<-1.0) return PI;
		else if(s>1.0) return 0;
		else return acos(s);
	}

	//float t = acos((a*a + b*b - c*c)/(2*a*b));
	//return t;
}


/*
It considers the orientation of (x1, y1) and (x2, y2) w.r.t. (x0,y0).
It returns negative if (x2, y2) is not clockwise orientation to (x1, y1).
It returns positive otherwise.
*/
float
GetVisualAngle2(float x1, float y1, float x2, float y2, 
			   float x0, float y0)
{
	int cw = ClockWise(x0, y0, x1, y1, x2, y2);
	if(cw>=0)
	{
		//return t;
		return GetVisualAngle(x1, y1, x2, y2, x0, y0);
	}
	else if(cw<0)
	{
		//return -t;
		return -GetVisualAngle(x1, y1, x2, y2, x0, y0);
	}
}

CParticleF
representativePosition(const vector<CParticleF>& cluster)
{
	CParticleF center(0, 0, 0);
	for(int i=0; i<cluster.size(); ++i)
	{
		center.m_X += cluster[i].m_X;
		center.m_Y += cluster[i].m_Y;
	}
	center.m_X /= (float)cluster.size();
	center.m_Y /= (float)cluster.size();

	CParticleF rep = cluster[0];
	float dmin = Distance(center, cluster[0]);
	for(int i=1; i<cluster.size(); ++i)
	{
		float d = Distance(center, cluster[i]);
		if(d < dmin)
		{
			dmin = d;
			rep = cluster[i];
		}
	}
	return rep;
}

float
areaTriangle(CParticleF& a, CParticleF& b, CParticleF& c)
{
	float A = Distance(a, b);
	float B = Distance(b, c);
	float C = Distance(c, a);
	float s = (A + B+ C)/2.0f;
	if(s<=A || s<=B || s<=C) return 0.0f;
	else
		return sqrt(s*(s-A)*(s-B)*(s-C));
}

/*
*/
bool
inside(const CParticleF& p, const vector<CParticleF>& pnts)
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
			(x < (float)(xj-xi) * (float)(y-yi) / (float)(yj-yi) + xi) )
			c = !c;
	}
	return c;
}

/*
compute the closest point from x on the line connecting p and q.
p and q have to be distinct.
*/
CParticleF Closest2Line(const CParticleF& p, const CParticleF& q, const CParticleF& x)
{
	float a = p.m_X;
	float b = q.m_X-p.m_X;
	float c = p.m_Y;
	float d = q.m_Y-p.m_Y;
	float t = (-a*b - c*d + b*x.m_X + d*x.m_Y)/(b*b + d*d);
	return CParticleF(a+b*t, c+d*t);
}

/*
compute the distance from x to the line connecting p and q.
p and q have to be distinct.
*/
float Distance2Line(const CParticleF& p, const CParticleF& q, const CParticleF& x)
{
	CParticleF y = Closest2Line(p, q, x);
	return Distance(x, y);
}

/*
compute the distance from x to the line segment connecting p and q.
If the closest point is outside the segment [p, q], then it selects the closer of the two end points.
p and q have to be distinct.
*/
float Distance2LineSegment(const CParticleF& p, const CParticleF& q, const CParticleF& x)
{
	float a = p.m_X;
	float b = q.m_X-p.m_X;
	float c = p.m_Y;
	float d = q.m_Y-p.m_Y;
	float t = (-a*b - c*d + b*x.m_X + d*x.m_Y)/(b*b + d*d);
	if(t>=0 && t<=1.0f)
	{
		CParticleF y(a+b*t, c+d*t);
		float d = Distance(x, y);
		return d;
	}
	else if(t>1.0f)
	{
		return Distance(q, x);
	}
	else
	{
		return Distance(p, x);
	}
}

float
polygonArea(vector<CParticleF>& vp)
{
	if(vp.size()<3) return 0.0f;

	float area = 0;
	for(int i=0, j=vp.size()-1; i<vp.size(); i++)
	{
		area += vp[j].m_X*vp[i].m_Y - vp[i].m_X*vp[j].m_Y;
		j=i;
	}
	return Abs(area) /2.0;
}

