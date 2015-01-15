#ifndef ___SZ_CONVOLUTION_TEMPLATE_H___
#define ___SZ_CONVOLUTION_TEMPLATE_H___

#include <vector>
using namespace std;

/*
So far, only 'extension by zero' has been implemented.
*/
enum ConvolutionBoundaryExtension {CBE_ZERO, CBE_COPY, CBE_MIRROR, CBE_WRAP};

template<class Item>
void
separableConvolution(vector<Item>& R,
					 const vector<Item>& A,
					 const vector<Item>& f,
					 ConvolutionBoundaryExtension extension,
					 int ndim,
					 const int* dims);

template<class Item>
void
separableVolumeConvolution(vector<Item>& R,
						   const vector<Item>& A,
						   const vector<Item>& f,
						   ConvolutionBoundaryExtension extension,
						   const int* dims);

template<class Item>
void
separableConvolution(vector<Item>& R,
					 const vector<Item>& A,
					 const vector<vector<Item>>& vf,
					 ConvolutionBoundaryExtension extension,
					 int ndim,
					 const int* dims);

template<class Item>
vector<Item>
makeGaussianFilter(Item mean, Item std, int length, bool bNormalize = true);

#include <../src/szConvolutionTemplate.cpp>

#endif /* ___SZ_CONVOLUTION_TEMPLATE_H___ */
