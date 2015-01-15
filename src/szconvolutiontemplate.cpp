#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>

template<class Item>
vector<Item>
makeGaussianFilter(Item mean, 
				   Item std, 
				   int length, 
				   bool bNormalize)
{
	vector<Item> filter(length);
	int offset = length/2;
	float sum = 0;
	int i;
	for(i=0; i<length; ++i)
	{
		Item val = exp(-((Item)(i-offset)*(i-offset)/(2*std*std)));
		sum += val;
		filter[i] = val;
	}
	if(bNormalize)
	{
		for(i=0; i<length; ++i)
		{
			filter[i] /= sum;
		}
	}
	return filter;
}

template<class Item>
void
OneDimensionConvolution(vector<Item>& B, 
						const vector<Item> &A, 
						const vector<Item>& f, 
						int n, 
						ConvolutionBoundaryExtension extension,
						int ndim, 
						const int* dims)
{
	int i, j;
	vector<int> vsub(ndim);
	int nvoxels = numberOfElements(ndim, dims);
	int offset = (f.size()-1)/2;
	int stride = 1;
	for(i=0; i<n; ++i)
	{
		stride *= dims[i];
	}

	Item val;
	for(i=0; i<nvoxels; ++i)
	{
		Item sum = 0;
		for(j=0; j<f.size(); ++j)
		{
			if(NeighborCheck(i, (j-offset)*stride, ndim, dims))
			{
				Item val = A[i + (j-offset)*stride];
				sum += val * f[j];
			}
		}
		B[i] = sum;
	}
}

template<class Item>
void
OneDimensionVolumeConvolution(vector<Item>& B, 
							  const vector<Item> &A, 
							  const vector<Item>& f, 
							  int direction, 
							  ConvolutionBoundaryExtension extension,
							  const int* dims)
{
	int i, j, k, m;
	int offset = (f.size()-1)/2;
	if(direction==0)
	{
		for(i=0; i<dims[2]; ++i) 
		{
			for(j=0; j<dims[1]; ++j) 
			{
				for(k=0; k<dims[0]; ++k) 
				{
					Item sum = 0;
					for(m=0; m<f.size(); ++m)
					{
						int x = k + m - offset;
						if(x>=0 && x<dims[0])
						{
							Item val = GetData3(A, x, j, i, dims[0], dims[1], dims[2], (Item)0);
							sum += val * f[m];
						}
					}
					SetData3(B, k, j, i, dims[0], dims[1], dims[2], sum);
				}
			}
		}
	}
	else if(direction==1)
	{
		for(i=0; i<dims[2]; ++i) 
		{
			for(j=0; j<dims[1]; ++j) 
			{
				for(k=0; k<dims[0]; ++k) 
				{
					Item sum = 0;
					for(m=0; m<f.size(); ++m)
					{
						int y = j + m - offset;
						if(y>=0 && y<dims[1])
						{
							Item val = GetData3(A, k, y, i, dims[0], dims[1], dims[2], (Item)0);
							sum += val * f[m];
						}
					}
					SetData3(B, k, j, i, dims[0], dims[1], dims[2], sum);
				}
			}
		}
	}
	else if(direction==2)
	{
		for(i=0; i<dims[2]; ++i) 
		{
			for(j=0; j<dims[1]; ++j) 
			{
				for(k=0; k<dims[0]; ++k) 
				{
					Item sum = 0;
					for(m=0; m<f.size(); ++m)
					{
						int z = i + m - offset;
						if(z>=0 && z<dims[2])
						{
							Item val = GetData3(A, k, j, z, dims[0], dims[1], dims[2], (Item)0);
							sum += val * f[m];
						}
					}
					SetData3(B, k, j, i, dims[0], dims[1], dims[2], sum);
				}
			}
		}
	}
}

template<class Item>
void
separableConvolution(vector<Item>& R,
					 const vector<Item>& A,
					 const vector<Item>& f,
					 ConvolutionBoundaryExtension extension,
					 int ndim,
					 const int* dims)
{
	vector<Item> T(A.size()); //temporary
	int n;
	for(n=0; n<ndim; ++n)
	{
		if(n == 0)
		{
			OneDimensionConvolution(R, A, f, n, extension, ndim, dims);
		}
		else if(n % 2)
		{
			OneDimensionConvolution(T, R, f, n, extension, ndim, dims);
		}
		else
		{
			OneDimensionConvolution(R, T, f, n, extension, ndim, dims);
		}
	}
	if(ndim % 2)
	{
		R = T;
	}
}

template<class Item>
void
separableConvolution(vector<Item>& R,
					 const vector<Item>& A,
					 const vector<vector<Item>>& vf,
					 ConvolutionBoundaryExtension extension,
					 int ndim,
					 const int* dims)
{

	vector<Item> T(A.size()); //temporary
	int n;
	for(n=0; n<ndim; ++n)
	{
		if(n == 0)
		{
			OneDimensionConvolution(R, A, vf[n], n, extension, ndim, dims);
		}
		else if(n % 2)
		{
			OneDimensionConvolution(T, R, vf[n], n, extension, ndim, dims);
		}
		else
		{
			OneDimensionConvolution(R, T, vf[n], n, extension, ndim, dims);
		}
	}
	if(ndim % 2)
	{
		R = T;
	}
}

template<class Item>
void
separableVolumeConvolution(vector<Item>& R,
						   const vector<Item>& A,
						   const vector<Item>& f,
						   ConvolutionBoundaryExtension extension,
						   const int* dims)
{
	vector<Item> T(A.size()); //temporary
	OneDimensionVolumeConvolution(R, A, f, 0, extension, dims);
	OneDimensionVolumeConvolution(T, R, f, 1, extension, dims);
	OneDimensionVolumeConvolution(R, T, f, 2, extension, dims);
}