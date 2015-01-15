#ifndef ___SZ_CPARTICLEF_H___
#define ___SZ_CPARTICLEF_H___
#include <cmath>
#include <limits>

class CParticleF
{
public:
	CParticleF(float x=0, float y=0, float z=0, float life=0)
	{
		m_X = x;
		m_Y = y;
		m_Z = z;
		m_Life = life;
	}
	bool operator == (const CParticleF& p) const
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
	bool operator < (const CParticleF& p) const
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
	bool operator <= (const CParticleF& p) const
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
	bool operator > (const CParticleF& p) const
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
	bool operator >= (const CParticleF& p) const
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
	float m_X;
	float m_Y;
	float m_Z;
	float m_Life;
};

inline CParticleF 
Difference(const CParticleF& p, const CParticleF& q)
{
	return CParticleF(q.m_X - p.m_X, q.m_Y - p.m_Y, q.m_Z - p.m_Z);
}

inline float
Length(const CParticleF& p)
{
	return sqrt(p.m_X*p.m_X + p.m_Y*p.m_Y + p.m_Z*p.m_Z);
}

inline CParticleF
UnitVector(const CParticleF& p)
{
	float len = Length(p);
	if(len <= std::numeric_limits<float>::epsilon() * 10.0f) len=1.0f; //too small...

	return CParticleF(p.m_X/len, p.m_Y/len, p.m_Z/len);
}

inline float
Distance(const CParticleF& p, const CParticleF& q)
{
	return sqrt((float)(p.m_X-q.m_X)*(p.m_X-q.m_X)+(p.m_Y-q.m_Y)*(p.m_Y-q.m_Y)+(p.m_Z-q.m_Z)*(p.m_Z-q.m_Z));
}

inline float
AngleXY(const CParticleF& p, const CParticleF& q)
{
	return atan2((float)(q.m_Y-p.m_Y), (float)(q.m_X-p.m_X));
}

class CParticleFL: public CParticleF
{
public:
	CParticleFL(float x=0, float y=0, float z=0, float life=0, int label=0, bool b=false): CParticleF(x, y, z, life)
	{
		m_Label = label;
		m_Flag = b;
	}
	CParticleFL(const CParticleF& p): CParticleF(p)
	{
		m_Label = 0;
		m_Flag = false;
	}
	int m_Label;
	bool m_Flag;
};

#endif /* ___SZ_CPARTICLE_H___ */