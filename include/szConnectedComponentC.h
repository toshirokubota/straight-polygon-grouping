#ifndef ___CONNECTED_COMPONENTC_H___
#define ___CONNECTED_COMPONENTC_H___

#include "szMexUtility.h"
#include "szMexUtilityTemplate.h"
#include "szMyNeighborOp.h"

/*
This is identical to ConnectedComponent.h except that:
the size of each component is returned in vcount.
*/

/*
CCA with general neighborhood
*/
template<class Item>
int
ConnectedComponentAnalysisEqual(vector<int>& C, 
                                const vector<Item>& V, 
                                const vector<int>& nbh,
								vector<int>& vcount,
                                Item value, 
                                int ndim, 
                                const int* dims);

/*
A wrapper for CCA
Use mode to specifiy a commonn neighborhood type.
*/
template<class Item>
int
ConnectedComponentAnalysisEqual(vector<int>& C, 
                                const vector<Item>& V, 
								vector<int>& vcount,
                                int mode,
                                Item value, 
                                int ndim, 
                                const int* dims);

/*
CCA with general neighborhood - size threshold is implemented
*/
template<class Item>
int
ConnectedComponentAnalysisEqualSizeThreshold(vector<int>& C, 
                                             const vector<Item>& V, 
                                             const vector<int>& nbh,
											 vector<int>& vcount,
                                             Item value, 
                                             int min_size,
                                             int max_size,
                                             int ndim, 
                                             const int* dims);

/*
A wrapper for CCA
Use mode to specifiy a commonn neighborhood type.
*/
template<class Item>
int
ConnectedComponentAnalysisEqualSizeThreshold(vector<int>& C, 
                                             const vector<Item>& V, 
											 vector<int>& vcount,
                                             int mode,
                                             Item value, 
                                             int min_size,
                                             int max_size,
                                             int ndim, 
                                             const int* dims);

/*
CCA with general neighborhood
*/
template<class Item>
int
ConnectedComponentAnalysisBigger(vector<int>& C, 
                                 const vector<Item>& V, 
                                 const vector<int>& nbh,
								 vector<int>& vcount,
                                 Item value, 
                                 int ndim, 
                                 const int* dims);

/*
A wrapper for CCA
Use mode to specifiy a commonn neighborhood type.
*/
template<class Item>
int
ConnectedComponentAnalysisBigger(vector<int>& C, 
                                 const vector<Item>& V, 
								 vector<int>& vcount,
                                 int mode,
                                 Item value, 
                                 int ndim, 
                                 const int* dims);


/*
CCA with general neighborhood
*/
template<class Item>
int
ConnectedComponentAnalysisSmaller(vector<int>& C, 
                                 const vector<Item>& V, 
                                 const vector<int>& nbh,
								 vector<int>& vcount,
                                 Item value, 
                                 int ndim, 
                                 const int* dims);

/*
A wrapper for CCA
Use mode to specifiy a commonn neighborhood type.
*/
template<class Item>
int
ConnectedComponentAnalysisSmaller(vector<int>& C, 
                                 const vector<Item>& V, 
								 vector<int>& vcount,
                                 int mode,
                                 Item value, 
                                 int ndim, 
                                 const int* dims);


#include "../src/szConnectedComponentC.cpp"

#endif /* ___CONNECTED_COMPONENTC_H___ */