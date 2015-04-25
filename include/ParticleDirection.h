#pragma once

struct ParticleDirection {
	float x;
	float y;
	float z;
	ParticleDirection(float x = 0, float y = 0, float z = 0) 
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	ParticleDirection(CParticleF& from, CParticleF& to)
	{
		x = to.m_X - from.m_X;
		y = to.m_Y - from.m_Y;
		z = to.m_Z - from.m_Z;
	}
};