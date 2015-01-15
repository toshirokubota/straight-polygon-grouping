#ifndef __SZ_CONTOUR_H__
#define __SZ_CONTOUR_H__
#include <mex.h>
#include <vector>
using namespace std;
#include <szParticleF.h>

struct Contour
{
	vector<CParticleF> points;
};

int
LoadContour(vector<Contour>& contours, const mxArray *prhs);

mxArray*
StoreContours(const vector<Contour>& contours);

#endif /* __SZ_CONTOUR_H__ */