#ifndef _SZ_MISC_OPERATIONS_H_
#define _SZ_MISC_OPERATIONS_H_
#include <vector>
using namespace std;
#include <szParticleF.h>
const double PI = 3.1415926536;

int GetOffsetX6(int n);

int GetOffsetY6(int n);

int GetOffsetZ6(int n);

int GetOffsetX10(int n);

int GetOffsetY10(int n);

int GetOffsetZ10(int n);

int GetOffsetX26(int n);

int GetOffsetY26(int n);

int GetOffsetZ26(int n);

vector<real> 
SymmetricExtension(const vector<real>& A, 
				   const int* dims);

void
CurvatureImage(vector<real>& K,
			   const vector<real>& U,
			   const int* dims);

float
Distance(const CParticleF& p1, const CParticleF& p2);

void
RemoveJunctions(vector<unsigned char>& I,
				const int* dims);

#endif /* _SZ_MISC_OPERATIONS_H_ */