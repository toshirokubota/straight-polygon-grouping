#ifndef ___LOCAL_EXTREMA_H___
#define ___LOCAL_EXTREMA_H___

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>

const int LocalExtremaStrict = 0;
const int LocalExtremaNonStrict = 1;

template<class Item>
void
LocalMaximum(vector<unsigned char>& M, 
             const vector<Item>& V, 
             const vector<unsigned char>& L,
             const vector<int>& nbh,
             bool strict,
             int ndim,
             const int* dims);

template<class Item>
void
LocalMinimum(vector<unsigned char>& M, 
             const vector<Item>& V, 
             const vector<unsigned char>& L,
             const vector<int>& nbh,
             bool strict,
             int ndim,
             const int* dims);

#include "../src/szLocalExtrema.cpp"

#endif /* ___LOCAL_EXTREMA_H___ */