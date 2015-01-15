#ifndef ___SZ_CPARTICLE_H___
#define ___SZ_CPARTICLE_H___
#include <cmath>

class CParticle
{
public:
	CParticle(int x=0, int y=0, int z=0, float life=0)
	{
		m_X = x;
		m_Y = y;
		m_Z = z;
		m_Life = life;
	}
	bool operator == (const CParticle& p) const
	{
		if(m_X == p.m_X && m_Y == p.m_Y && m_Z == p.m_Z)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator < (const CParticle& p) const
	{
		if(this->m_Life < p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator <= (const CParticle& p) const
	{
		if(this->m_Life <= p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator > (const CParticle& p) const
	{
		if(this->m_Life > p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator >= (const CParticle& p) const
	{
		if(this->m_Life >= p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
public:
	int m_X;
	int m_Y;
	int m_Z;
	float m_Life;
};

inline CParticle
Difference(const CParticle& p, const CParticle& q)
{
	return CParticle(q.m_X - p.m_X, q.m_Y - p.m_Y, q.m_Z - p.m_Z);
}


inline float
Distance(const CParticle& p, const CParticle& q)
{
	return sqrt((float)(p.m_X-q.m_X)*(p.m_X-q.m_X)+(p.m_Y-q.m_Y)*(p.m_Y-q.m_Y)+(p.m_Z-q.m_Z)*(p.m_Z-q.m_Z));
}

inline float
AngleXY(const CParticle& p, const CParticle& q)
{
	return atan2((float)(q.m_Y-p.m_Y), (float)(q.m_X-p.m_X));
}


#endif /* ___SZ_CPARTICLE_H___ */