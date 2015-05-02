#pragma once;
#include <MovingParticle.h>
#include <mex.h>
#include <ParticleSimulator.h>

class ParticleSimulatorV2: public ParticleSimulator
{
public:
	virtual bool Simulate(float endtime, float delta, bool bdebug);
private:
};

