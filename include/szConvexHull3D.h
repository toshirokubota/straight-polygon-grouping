#ifndef _SZ_CONVEX_HULL_3D_H_
#define _SZ_CONVEX_HULL_3D_H_

#include <vector>
using namespace std;

void
preConvexHull(vector<unsigned char>& S,   //OUTPUT - must be initialized to zero
              const vector<unsigned char>& C,   //INPUT - segmentation map
              const int* dims);

void
postConvexHull(vector<unsigned char>& B,    //OUTPUT - must be initialized to zero
               const vector<unsigned char>& C,    //INPUT - convex-hull map
               const vector<unsigned char>& L,   //INPUT - foreground segmentation map
               const int* dims);             //volume size


bool
ConvexHull3D (vector<unsigned char>& S, 
              const vector<unsigned char>& L, 
              int ndim, 
              const int* dims);

#endif /* _SZ_CONVEX_HULL_3D_H_ */