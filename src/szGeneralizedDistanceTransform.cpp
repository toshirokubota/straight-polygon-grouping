#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szDistanceTransform.h>
#include <limits>
using namespace std;

void
GeneralizedDistanceTransform(vector<double>& D, 
                        const vector<double>& L, 
						double maxValue, //maximum valued allowed in L
                        int ndim, 
                        const int* dims) 
{
	int i,n;

	for(i=0; i<D.size(); ++i) 
	{
		if(L[i]<=0)
			D[i]=std::numeric_limits<double>::infinity();
		else
			D[i]=1.0/L[i] - 1.0/maxValue;
	}

	int nvoxels=numberOfElements(ndim,dims);
	int stride=1;
	vector<int> vsub(ndim); //subscript buffer
	for(n=0; n<ndim; ++n) 
	{
		vector<double> buffer(dims[n]);
		int offset=0;
		for(int m=0; m<nvoxels/dims[n]; m++) 
		{
			//compute the offset
			int m2=m;
			int k;
			for(k=0; k<ndim; ++k) 
			{
				if(k!=n) 
				{
					vsub[k]=m2 % dims[k];
					m2/=dims[k];
				}
				else
					vsub[k]=0;
			}
			int offset=0;
			int stride2=1;
			for(k=0; k<ndim; ++k) 
			{
				offset+=vsub[k]*stride2;
				stride2*=dims[k];
			}

			//copy relevant line of data
			for(k=0; k<dims[n]; ++k)
				buffer[k]=GetData(D,offset+k*stride,(double)0);

			//do the computation
			vector<double> dst = dt_Felzenszwalb(buffer,dims[n]);

			//update the distance
			for(k=0; k<dims[n]; ++k)
				SetData(D,offset+k*stride,dst[k]);
		}
		stride*=dims[n];
	}

	//take the square root of the distance square
	for(n=0; n<D.size(); ++n)
		D[n]=sqrt(D[n]); //no array limit checking...
}

void
GeneralizedDistanceTransformF(vector<float>& D, 
                        const vector<float>& L, 
						float maxValue, //maximum valued allowed in L
                        int ndim, 
                        const int* dims) 
{
	int i,n;

	for(i=0; i<D.size(); ++i) 
	{
		if(L[i]<=0)
			D[i]=std::numeric_limits<float>::infinity();
		else
			D[i]=1.0/L[i] - 1.0/maxValue;
	}

	int nvoxels=numberOfElements(ndim,dims);
	int stride=1;
	vector<int> vsub(ndim); //subscript buffer
	for(n=0; n<ndim; ++n) 
	{
		vector<double> buffer(dims[n]);
		int offset=0;
		for(int m=0; m<nvoxels/dims[n]; m++) 
		{
			//compute the offset
			int m2=m;
			int k;
			for(k=0; k<ndim; ++k) 
			{
				if(k!=n) 
				{
					vsub[k]=m2 % dims[k];
					m2/=dims[k];
				}
				else
					vsub[k]=0;
			}
			int offset=0;
			int stride2=1;
			for(k=0; k<ndim; ++k) 
			{
				offset+=vsub[k]*stride2;
				stride2*=dims[k];
			}

			//copy relevant line of data
			for(k=0; k<dims[n]; ++k)
				buffer[k]=(double)GetData(D,offset+k*stride,(float)0);

			//do the computation
			vector<double> dst = dt_Felzenszwalb(buffer,dims[n]);

			//update the distance
			for(k=0; k<dims[n]; ++k)
				SetData(D,offset+k*stride,(float)dst[k]);
		}
		stride*=dims[n];
	}

	//take the square root of the distance square
	for(n=0; n<D.size(); ++n)
		D[n]=sqrt(D[n]); //no array limit checking...
}
