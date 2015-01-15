#ifndef ___MY_MEX_FILE_IO_H____
#define ___MY_MEX_FILE_IO_H____
#include <mex.h>
#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

#include <szMexUtility.h>

template<class Item>
int
LoadData(vector<Item>& A, const mxArray *prhs, mxClassID& class_id, int& ndim, const int** dims);

template<class Item>
int
ReadScalar(Item& A, const mxArray *prhs, mxClassID& class_id);

template<class Item>
mxArray*
StoreData(vector<Item>& vdata, mxClassID class_id, int ndims, const int* dims);

template<class Item>
mxArray*
StoreScalar(Item data, mxClassID class_id);

template<class Item>
mxArray*
StoreDataCell(vector<vector<Item>>& vdata, mxClassID class_id, int ndim, const int* dims, int ncols=1);

#include "../src/MexFileIO.cpp"

#endif /* ___MY_MEX_FILE_IO_H____ */