#ifndef ___MY_MEX_UTILITY_H____
#define ___MY_MEX_UTILITY_H____

#include <string.h>
#include <math.h>
#include <vector>
using namespace std;
#include <myDataType.h>
//typedef const int mwSize;

#define Max(a,b) ((a)>(b) ? (a):(b))
#define Min(a,b) ((a)<(b) ? (a):(b))
#define Abs(a) ((a)>0 ? (a):(-(a)))
#define Round(a) ((int)((a)+.5))
const double INFINITY = 1.0E300;  //for a practical purpose, this is enough.
const float INFINITYF = (float)1.0E38; //for a practical purpose, this is enough.
const int MAX_INTEGER = 0x7FFFFFFF;
const int MIN_INTEGER = 0x80000000;
const int MAX_UNSIGNED_INTEGER = 0xFFFFFFFF;
const int MAX_SHORT = 0x7FFF;
const int MIN_SHORT = 0x8000;
const int MAX_UNSIGNED_SHORT = 0xFFFF;

const int MyMexAllocationError = -1;
const int MyMexUnsupportedClassError = -2;
const int MyMexNotEnoughDataError = -3;

int 
Mod(int a, int b);

int
numberOfElements(int ndim, const int* dims);

int
Sub2Ind(int x, int y, int z, const int* dims);

void
Ind2Sub(int& x, int& y, int& z, int ind, const int* dims);

int
Sub2Ind(const vector<int>& vsub, int ndim, const int* dims);

vector<int>
Ind2Sub(int ind, int ndim, const int* dims);

void
Ind2Sub(vector<int>& vsub, int ind, int ndim, const int* dims);

int
Sub2Ind(const vector<int>& vsub, int ndim, int* dims);

vector<int>
Ind2Sub(int ind, int ndim, int* dims);

int
Sub2IndCentered(const vector<int>& vsub, int ndim, const int* dims);

vector<int>
Ind2SubCentered(int ind, int ndim, const int* dims);

int
Sub2IndCentered(const vector<int>& vsub, int ndim, int* dims);

vector<int>
Ind2SubCentered(int ind, int ndim, int* dims);

void
getNewDimension(int* newdims,
				vector<float>& vsize_new,
				vector<float>& vsize_old,
				int ndim,
				const int* olddims);

#endif
