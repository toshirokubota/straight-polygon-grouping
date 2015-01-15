#ifndef ___MY_NEIGHBOR_OP_H____
#define ___MY_NEIGHBOR_OP_H____

#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

#include "szMexUtility.h"
#include "szMexUtilityTemplate.h"

//neighborhood type
const int NeighborhoodFour = 1;
const int NeighborhoodEight = 2;

const int XOffset[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1};
const int YOffset[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
const int ZOffset[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int NumNeighbors = sizeof(XOffset)/sizeof(XOffset[0]);

const int XOffset6[] = {0, -1, 1, 0, 0, 0};
const int YOffset6[] = {-1, 0, 0, 1, 0, 0};
const int ZOffset6[] = {0, 0, 0, 0, -1, 1};
const int NumNeighbors6 = sizeof(XOffset6)/sizeof(XOffset6[0]);

const int XOffset10[] = {0, -1, 1, 0, 0, 0, 1, -1, 1, -1};
const int YOffset10[] = {-1, 0, 0, 1, 0, 0, 1, 1, -1, -1};
const int ZOffset10[] = {0, 0, 0, 0, -1, 1, 0, 0, 0, 0};
const int NumNeighbors10 = sizeof(XOffset10)/sizeof(XOffset10[0]);

const int XOffset18[] = {-1, 1, 0, 0, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1, 0, 0, 0, 0};
const int YOffset18[] = {0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0, -1, -1, 1, 1};
const int ZOffset18[] = {0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, -1, 1, -1, 1, -1, 1};
const int NumNeighbors18 = sizeof(XOffset18)/sizeof(XOffset18[0]);

const int XOffset8[] = {-1, 0, 1, -1, 1, -1, 0, 1};
const int YOffset8[] = {-1, -1, -1, 0, 0, 1, 1, 1};
const int NumNeighbors8 = sizeof(XOffset8)/sizeof(XOffset8[0]);

const int XOffset4[] = {0, -1, 1, 0};
const int YOffset4[] = {-1, 0, 0, 1};
const int NumNeighbors4 = sizeof(XOffset4)/sizeof(XOffset4[0]);

bool
BoundaryCheck(int index, 
              int ndim, 
              const int* dims);

bool
NeighborCheck(int p,
              int off,
              int ndim, 
              const int* dims);

bool
NeighborCheck(const int* p,
              const int* q,
              int ndim, 
              const int* dims);

bool
NeighborCheck(vector<int>::const_iterator p,
              vector<int>::const_iterator q,
              int ndim, 
              const int* dims);

vector<int>
makeNeighborhood10(int ndim, const int* dims);

vector<int>
MakeFourNeighborhood(int ndim,
                     const int* dims);

vector<int>
MakeCausalFourNeighborhood(int ndim,
                     const int* dims);

vector<int>
MakeEightNeighborhood(int ndim,
                     const int* dims);

vector<int>
MakeNineNeighborhood(int ndim,
                     const int* dims);

vector<int>
MakeAntiCausalNeighborhood(int ndim,
                           const int* dims);

vector<int>
MakeCausalNeighborhood(int ndim,
                           const int* dims);

vector<int>
MakeNeighborhood(const int* width,
				 int ndim,
				 const int* dims);

vector<int>
MakeNeighborhood(int width,
				 int ndim,
				 const int* dims);

#endif /* ___MY_NEIGHBOR_OP_H____ */