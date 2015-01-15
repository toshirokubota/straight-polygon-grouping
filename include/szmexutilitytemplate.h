#ifndef ___MY_MEX_UTILITY_TEMPLATE_H____
#define ___MY_MEX_UTILITY_TEMPLATE_H____

#include <fstream>
#include <string>
#include <cmath>
#include <vector>
using namespace std;

#include "szMexUtility.h"
#include "szParticleF.h"

//Generic data assess
template<class Item>
Item
GetData(const vector<Item>& A, int i, bool& success);

template<class Item>
Item
GetData(const vector<Item>& A, int i, Item defval);

template<class Item>
bool
SetData(vector<Item>& A, int i, const Item val);

//3D data access
template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, bool& success);

template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, Item defval);

template<class Item>
bool
SetData3(vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, const Item val);

//2D data access
template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, bool& success);

template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, Item defval);

template<class Item>
bool
SetData2(vector<Item>& A, int x, int y, int xD, int yD, const Item val);

//N-D data access
template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, Item& defval);

template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, bool& success);

template<class Item>
bool
SetDataN(vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, const Item& val);

template<class T>
bool
SetVoxel(vector<T>& A, 
		 const CParticleF& p,
		 T value,
		 const int* dims);

template<class T>
T
GetVoxel(const vector<T>& A, 
		 const CParticleF& p,
		 T defaultValue,
		 const int* dims);


#include "../src/szMexUtilityTemplate.cpp"

#endif /* ___MY_MEX_UTILITY_TEMPLATE_H____ */