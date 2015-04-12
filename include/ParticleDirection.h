#pragma once

struct ParticleDirection {
	float x;
	float y;
	float z;
	ParticleDirection(float x = 0, float y = 0, float z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
};