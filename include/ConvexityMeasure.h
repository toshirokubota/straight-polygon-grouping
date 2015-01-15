#ifndef ___CONVEXITY_MEASURE_H___
#define ___CONVEXITY_MEASURE_H___

#include <szParticleF.h>
#include <vector>
#include <Triangulation.h>

float
ConvexityMeasure(const std::vector<CParticleF>& points, float maxGap);

float
ConvexityMeasure(Triangulation::Triangulator& trmap);

#endif /* ___CONVEXITY_MEASURE_H___ */