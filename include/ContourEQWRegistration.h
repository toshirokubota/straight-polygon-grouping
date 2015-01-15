#ifndef __CONTOUR_EQW_REGISTRATION_H___
#define __CONTOUR_EQW_REGISTRATION_H___
#include <mex.h>
#include <vector>
using namespace std;
#include <ContourEQW.h>
#include <myDataType.h>

ContourEQW Interpolate(const ContourEQW& a);

CParticleF Centroid(const ContourEQW& c);

/*
Using a rigid transformation, register contour B to contour A.
The resulting contour has the same length with B. 
Only X and Y are set and coefficients are set to 0.
*/
ContourEQW
AffineRegistration(const ContourEQW& a, const ContourEQW& b);

ContourEQW
AffineRegistration(const ContourEQW& a, const ContourEQW& b, vector<real>& params);

/*
Using a rigid transformation, register contour B to contour A.
The resulting contour has the same length with B. No translation is involved.
Only X and Y are set and coefficients are set to 0.
*/
ContourEQW
EuclideanRegistration(const ContourEQW& a, const ContourEQW& b);

ContourEQW
EuclideanRegistration(const ContourEQW& a, const ContourEQW& b, vector<real>& params);

/*
Using a Procrustes transformation, register contour B to contour A.
Transform consists of translation, uniform scaling, and rotation.
The resulting contour has the same length with B. 
Only X and Y are set and coefficients are set to 0.
*/
ContourEQW
ProcrustesRegistration(const ContourEQW& a, const ContourEQW& b, vector<real>& z);

ContourEQW
ProcrustesRegistration(const ContourEQW& a, const ContourEQW& b);

#endif /*  __CONTOUR_EQW_H___ */