#include <algorithm>
using namespace std;

//Generic data access
template<class Item>
Item
GetData(const vector<Item>& A, int i, bool& success) {
	if (i>=0 && i<A.size()) {
		success=true;
		return A[i];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData(const vector<Item>& A, int i, Item defval) {
	if (i>=0 && i<A.size()) {
		return A[i];
	}
	else {
		return (Item)defval;
	}
}

template<class Item>
bool
SetData(vector<Item>& A, int i, const Item val) {
	if (i>=0 && i<A.size()) {
		A[i]=val;
		return true;
	}
	else {
		return false;
	}
}

//3D data access
template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, bool& success) {
	if (x>=0 && x<xD && y>=0 && y<yD && z>=0 && z<zD) {
		success=true;
		return A[z*xD*yD+y*xD+x];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData3(const vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, Item defval) {
	if (x>=0 && x<xD && y>=0 && y<yD && z>=0 && z<zD) {
		return A[z*xD*yD+y*xD+x];
	}
	else {
		return defval;
	}
}

template<class Item>
bool
SetData3(vector<Item>& A, int x, int y, int z, int xD, int yD, int zD, Item val) {
	if (x>=0 && x<xD && y>=0 && y<yD && z>=0 && z<zD) {
		A[z*xD*yD+y*xD+x]=val;
		return true;
	}
	else {
		return false;
	}
}

//2D data access
template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, bool& success) {
	if (x>=0 && x<xD && y>=0 && y<yD) {
		success=true;
		return A[y*xD+x];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
Item
GetData2(const vector<Item>& A, int x, int y, int xD, int yD, Item defval) {
	if (x>=0 && x<xD && y>=0 && y<yD) {
		return A[y*xD+x];
	}
	else {
		return defval;
	}
}

template<class Item>
bool
SetData2(vector<Item>& A, int x, int y, int xD, int yD, Item val) {
	if (x>=0 && x<xD && y>=0 && y<yD) {
		A[y*xD+x]=val;
		return true;
	}
	else {
		return false;
	}
}

//N-D data access
template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, Item defval) {
	int id=0;
	int skip=1;
	bool within=true;
	for(int i=0; i<ndim; ++i) {
		id+=vsub[i]*skip;
		skip*=dims[i];
		if(vsub[i]<0 || vsub[i]>=dims[i]) {
			within=false;
			break;
		}
	}

	if (within && id>=0 && id<skip) {
		return A[id];
	}
	else {
		return defval;
	}
}

template<class Item>
Item
GetDataN(const vector<Item>& A, const vector<int> vsub, const int* dims, int ndim, bool& success) {
	int id=0;
	int skip=1;
	bool within=true;
	for(int i=0; i<ndim; ++i) {
		id+=vsub[i]*skip;
		skip*=dims[i];
		if(vsub[i]<0 || vsub[i]>=dims[i]) {
			within=false;
			break;
		}
	}

	if (within && id>=0 && id<skip) {
		success=true;
		return A[id];
	}
	else {
		success=false;
		return (Item)0;
	}
}

template<class Item>
bool
SetDataN(vector<Item>& A, vector<int> vsub, const int* dims, int ndim, Item val) 
{
	int id=0;
	int skip=1;
	bool within=true;
	for(int i=0; i<ndim; ++i) 
	{
		id+=vsub[i]*skip;
		skip*=dims[i];
		if(vsub[i]<0 || vsub[i]>=dims[i]) 
		{
			within=false;
			break;
		}
	}

	if (within && id>=0 && id<skip) 
	{
		A[id]=val;
		return true;
	}
	else 
	{
		return false;
	}
}

template<class T>
bool
SetVoxel(vector<T>& A, 
		 const CParticleF& p,
		 T value,
		 const int* dims)
{
	return SetData3(A, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], value);
}

template<class T>
T
GetVoxel(const vector<T>& A, 
		 const CParticleF& p,
		 T defaultValue,
		 const int* dims)
{
	return GetData3(A, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], defaultValue);
}

