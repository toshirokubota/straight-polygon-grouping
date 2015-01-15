
#include <szDistanceTransform.h>
#include <szMyNeighborOp.h>
#define INFINITY_FELZENSWALB 1.0E300

vector<double>
makeQuasiEuclideanDistanceWeight(const vector<int>& nbh, 
                                 int ndim, 
                                 const int* dims) {
  vector<double> weight;
  for(int k=0; k<nbh.size(); ++k) {
    int ind=nbh[k];
    vector<int> vsub = Ind2SubCentered(ind,ndim,dims);
    double sum=0;
    for(int i=0; i<ndim; ++i) {
      sum+=vsub[i]*vsub[i];
    }
    weight.push_back(sqrt(sum));
  }

  for(int i=0; i<weight.size(); ++i) {
    //mexPrintf("makeQuasiEuclideanDistanceWeight: %f\n", weight[i]);
  }

  return weight;
}

vector<double>
makeCityBlockDistanceWeight(const vector<int>& nbh, 
                            int ndim, const 
                            int* dims) {
  vector<double> weight;
  for(int k=0; k<nbh.size(); ++k) {
    int ind=nbh[k];
    vector<int> vsub = Ind2SubCentered(ind,ndim,dims);
    double sum=0;
    for(int i=0; i<ndim; ++i) {
      sum+=Abs(vsub[i]);
    }
    weight.push_back(sum);
  }

  for(int i=0; i<weight.size(); ++i) {
    //mexPrintf("makeQuasiEuclideanDistanceWeight: %f\n", weight[i]);
  }

  return weight;
}

vector<double>
makeCheckerBoardDistanceWeight(const vector<int>& nbh, 
                               int ndim, 
                               const int* dims) {
  vector<double> weight;
  for(int k=0; k<nbh.size(); ++k) {
    int ind=nbh[k];
    vector<int> vsub = Ind2SubCentered(ind,ndim,dims);
    double maxv=0;
    for(int i=0; i<ndim; ++i) {
      maxv=Max(Abs(vsub[i]),maxv);
    }
    weight.push_back(maxv);
  }

  for(int i=0; i<weight.size(); ++i) {
    //mexPrintf("makeQuasiEuclideanDistanceWeight: %f\n", weight[i]);
  }

  return weight;
}

void
DistanceTransformNonEuclid(vector<double>& D, 
                           const vector<unsigned char>& L, 
                           const vector<int>& leadingNbh,
                           const vector<int>& trailingNbh,
                           const vector<double>& leadingWgt,
                           const vector<double>& trailingWgt,
                           int ndim, 
                           const int* dims) {

  int i,j;
  int nvoxels=numberOfElements(ndim,dims);

  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoordsL;
  for(i=0; i<leadingNbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(leadingNbh[i],ndim,dims);
    vcoordsL.insert(vcoordsL.end(),vsub.begin(),vsub.end());
  }
  vector<int> vcoordsT;
  for(i=0; i<leadingNbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(trailingNbh[i],ndim,dims);
    vcoordsT.insert(vcoordsT.end(),vsub.begin(),vsub.end());
  }

  //first pass - from upper left to lower right
  for(i=0; i<nvoxels; ++i) {
    double v1=GetData(D,i,(double)0);
    double minV=INFINITY;
    vector<int> vsub = Ind2Sub(i,ndim,dims);
    for(j=0; j<leadingNbh.size(); ++j) {
      if(NeighborCheck(vsub.begin(),vcoordsL.begin()+j*ndim,ndim,dims)) {
        int k=i+leadingNbh[j];
        double v2=D[k]; //trust BoundaryCheck - no array boundary check here
        minV=Min(minV,v2+leadingWgt[j]);
      }
    }
    if(v1>minV)
      SetData(D,i,minV);
  }

  //second pass - from lower right to upper left
  for(i=nvoxels-1; i>=0; --i) {
    double v1=GetData(D,i,(double)0);
    double minV=INFINITY;
    vector<int> vsub = Ind2Sub(i,ndim,dims);
    for(j=0; j<trailingNbh.size(); ++j) {
      if(NeighborCheck(vsub.begin(),vcoordsT.begin()+j*ndim,ndim,dims)) {
        int k=i+trailingNbh[j];
        double v2=D[k]; //trust BoundaryCheck - no array boundary check here
        minV=Min(minV,v2+trailingWgt[j]);
     }
    }
    if(v1>minV)
      SetData(D,i,minV);
  }
}

/* dt of 1d function using squared distance 
The original implementation is credited to

Distance Transforms of Sampled Functions
Pedro F. Felzenszwalb and Daniel P. Huttenlocher
Cornell Computing and Information Science TR2004-1963

*/

#define dt_square(a) ((a)*(a))
vector<double>  
dt_Felzenszwalb(const vector<double>& f, int n) 
{
  vector<double> d(n);
  vector<int> v(n);
  vector<double> z(n+1);
  int k = 0;
  v[0] = 0;
  z[0] = -INFINITY_FELZENSWALB;
  z[1] = +INFINITY_FELZENSWALB;
  int q;
  for (q = 1; q <= n-1; q++) {
    double s  = ((f[q]+dt_square(q))-(f[v[k]]+dt_square(v[k])))/(2*q-2*v[k]);
    while (s <= z[k]) {
      k--;
      s  = ((f[q]+dt_square(q))-(f[v[k]]+dt_square(v[k])))/(2*q-2*v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k+1] = +INFINITY_FELZENSWALB;
  }
  
  k = 0;
  for (q = 0; q <= n-1; q++) {
    while (z[k+1] < q)
      k++;
    d[q] = dt_square(q-v[k]) + f[v[k]];
  }
  return d;
}

void
DistanceTransformEuclid(vector<double>& D, 
                        const vector<unsigned char>& L, 
                        int ndim, 
                        const int* dims) 
{
	int i,n;

	for(i=0; i<D.size(); ++i) 
	{
		if(L[i])
			D[i]=0;
		else
			D[i]=INFINITY_FELZENSWALB;
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
DistanceTransform(vector<double>& D, 
                  const vector<unsigned char>& L, 
                  int mode,
                  int ndim, 
                  const int* dims) {
  int i;
  int nvoxels=numberOfElements(ndim,dims);
  for(i=0; i<nvoxels; ++i) {
    if(GetData(L,i,(unsigned char)0)) {
      SetData(D,i,(double)0);
    }
    else {
      SetData(D,i,INFINITY_FELZENSWALB);
    }
  }

  vector<int> leadingNbh=MakeCausalNeighborhood(ndim,dims);
  vector<int> trailingNbh=MakeAntiCausalNeighborhood(ndim,dims);;
  vector<double> leadingWgt;
  vector<double> trailingWgt;

  switch(mode) {
  case DistanceTransformMode_Euclid:
    DistanceTransformEuclid(D,L,ndim,dims);
    break;
  case DistanceTransformMode_QuasiEuclid:
    leadingWgt = makeQuasiEuclideanDistanceWeight(leadingNbh,ndim,dims);
    trailingWgt = makeQuasiEuclideanDistanceWeight(trailingNbh,ndim,dims);
    DistanceTransformNonEuclid(D,L,
      leadingNbh,trailingNbh,
      leadingWgt,trailingWgt,
      ndim,dims);
    break;
  case DistanceTransformMode_CityBlockEuclid:
    leadingWgt = makeCityBlockDistanceWeight(leadingNbh,ndim,dims);
    trailingWgt = makeCityBlockDistanceWeight(trailingNbh,ndim,dims);
    DistanceTransformNonEuclid(D,L,
      leadingNbh,trailingNbh,
      leadingWgt,trailingWgt,
      ndim,dims);
    break;
  case DistanceTransformMode_CheckerBoardEuclid:
    leadingWgt = makeCheckerBoardDistanceWeight(leadingNbh,ndim,dims);
    trailingWgt = makeCheckerBoardDistanceWeight(trailingNbh,ndim,dims);
    DistanceTransformNonEuclid(D,L,
      leadingNbh,trailingNbh,
      leadingWgt,trailingWgt,
      ndim,dims);
    break;
  }
}

void
DistanceTransformF(vector<float>& Df, 
                  const vector<unsigned char>& L, 
                  int mode,
                  int ndim, 
                  const int* dims) {
  int i;
  int nvoxels=numberOfElements(ndim,dims);
  vector<double> D(nvoxels);
  for(i=0; i<nvoxels; ++i) {
    if(GetData(L,i,(unsigned char)0)) {
      SetData(D,i,(double)0);
    }
    else {
      SetData(D,i,INFINITY_FELZENSWALB);
    }
  }

  vector<int> leadingNbh=MakeCausalNeighborhood(ndim,dims);
  vector<int> trailingNbh=MakeAntiCausalNeighborhood(ndim,dims);;
  vector<double> leadingWgt;
  vector<double> trailingWgt;

  switch(mode) {
  case DistanceTransformMode_Euclid:
    DistanceTransformEuclid(D,L,ndim,dims);
    break;
  case DistanceTransformMode_QuasiEuclid:
    leadingWgt = makeQuasiEuclideanDistanceWeight(leadingNbh,ndim,dims);
    trailingWgt = makeQuasiEuclideanDistanceWeight(trailingNbh,ndim,dims);
    DistanceTransformNonEuclid(D,L,
      leadingNbh,trailingNbh,
      leadingWgt,trailingWgt,
      ndim,dims);
    break;
  case DistanceTransformMode_CityBlockEuclid:
    leadingWgt = makeCityBlockDistanceWeight(leadingNbh,ndim,dims);
    trailingWgt = makeCityBlockDistanceWeight(trailingNbh,ndim,dims);
    DistanceTransformNonEuclid(D,L,
      leadingNbh,trailingNbh,
      leadingWgt,trailingWgt,
      ndim,dims);
    break;
  case DistanceTransformMode_CheckerBoardEuclid:
    leadingWgt = makeCheckerBoardDistanceWeight(leadingNbh,ndim,dims);
    trailingWgt = makeCheckerBoardDistanceWeight(trailingNbh,ndim,dims);
    DistanceTransformNonEuclid(D,L,
      leadingNbh,trailingNbh,
      leadingWgt,trailingWgt,
      ndim,dims);
    break;
  }

  for(i=0; i<nvoxels; ++i) {
    Df[i]=(float)D[i];
  }
}

void
DistanceTransformInv(vector<double>& D, 
					  const vector<unsigned char>& L, 
					  int mode,
					  int ndim, 
					  const int* dims) {

  vector<unsigned char> invL(L.size());
  int i;
  for(i=0; i<L.size(); ++i) {
    if(L[i]==0)
      invL[i]=1;
    else
      invL[i]=0;
  }

  DistanceTransform(D,invL,mode,ndim,dims);
}

void
DistanceTransformFInv(vector<float>& D, 
					  const vector<unsigned char>& L, 
					  int mode,
					  int ndim, 
					  const int* dims) {

  vector<unsigned char> invL(L.size());
  int i;
  for(i=0; i<L.size(); ++i) {
    if(L[i]==0)
      invL[i]=1;
    else
      invL[i]=0;
  }

  DistanceTransformF(D,invL,mode,ndim,dims);
}
