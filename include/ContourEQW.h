#ifndef __CONTOUR_EQW_H___
#define __CONTOUR_EQW_H___
#include <mex.h>
#include <vector>
using namespace std;
#include <szContour.h>
#include <myDataType.h>

struct ContourEQW
{
public:
	ContourEQW(int n=0);
	ContourEQW(Contour c);
	ContourEQW operator =(Contour c);
	void Smoothen(int niter, float lambda=1.0, float tau=1.0, bool bClosed=true, bool bFixedEnds=false);
	vector<float> Curvatures() const;
	float Orientation(int i) const;
	float Curvature(int i) const;
	float EvaluateX(int i, float t) const;
	float EvaluateY(int i, float t) const;
	ContourEQW extract(int i, int j) const;
	ContourEQW reverse() const;
	int size() const {return (int)X.size();}

	vector<float> X;
	vector<float> Y;
	vector<float> A;
	vector<float> B;
	vector<float> C;
	vector<float> D;
	vector<float> Strength;
	float G[3][3];
	float H[3][3];

	vector<float> _FindOptimumPoint(int n0, int n2) const;
};

int
LoadContourEQW(vector<ContourEQW>& contours, const mxArray *prhs);

mxArray*
StoreContoursEQW(const vector<ContourEQW>& contours);

/*
Extract a portion of a contour (first - last, inclusive) and returns it
as another contour. All params within the segment are kept intact.
*/
ContourEQW
extractContourFragment(ContourEQW& c, int first, int last);

#endif /*  __CONTOUR_EQW_H___ */