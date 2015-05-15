#pragma once
#include <ParticleDirection.h>

class MovingFront
{
public:
	MovingFront()
	{
		speed = 0.0f;
	}
	MovingFront(ParticleDirection dir, float speed = 1.0f)
	{
		this->dir = dir;
		this->speed = speed;
	}
	ParticleDirection normal() const {
		return ParticleDirection(-dir.y, dir.x);
	}
	ParticleDirection dir;
	float speed;
};
