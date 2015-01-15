#ifndef ___CIRCLE_FITTING_H___
#define ___CIRCLE_FITTING_H___
#include <vector>
#include <szParticleF.h>

bool fitCircle(const std::vector<CParticleF>& points,
				double& radius,
				double& centerx, 
				double & centery,
				int maxIter=10,
				double rate = 0.1);

double fittingCircleError(const std::vector<CParticleF>& points,
				double radius,
				double centerx, 
				double centery);

#endif /* ___CIRCLE_FITTING_H___ */