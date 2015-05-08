#pragma once
#include <ParticleDirection.h>

class MovingFront
{
public:
	MovingFront(ParticleDirection dir, float speed = 1.0f)
	{
		this->dir = dir;
		this->speed = speed;
	}
	ParticleDirection dir;
	float speed;
};
