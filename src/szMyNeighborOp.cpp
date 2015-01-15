#include <szMyNeighborOp.h>
#ifdef MEX_DLL
#include <mex.h>
#endif

bool
BoundaryCheck(const vector<int> vsub, int ndim, const int* dims) {
  for(int i=0; i<ndim; ++i) {
    if(vsub[i]<0 || vsub[i]>=dims[i]) {
      return false;
    }
  }
  return true;
}

bool
NeighborCheck(int p,
              int offset,
              int ndim, 
              const int* dims) {
  vector<int> vsub = Ind2Sub(p, ndim, dims);
  vector<int> voff = Ind2SubCentered(offset, ndim, dims);
  if(voff.size()>ndim)
    return false; //index overflow

  for(int i=0; i<ndim; ++i) {
    if(vsub[i]+voff[i]<0 || vsub[i]+voff[i]>=dims[i]) {
      return false;
    }
  }
  return true;
}

bool
NeighborCheck(const int* p,
              const int* q,
              int ndim, 
              const int* dims) {
  for(int i=0; i<ndim; ++i) {
    int k=p[i]+q[i];
    if(k<0 || k>=dims[i]) {
      return false;
    }
  }
  return true;
}

bool
NeighborCheck(vector<int>::const_iterator p,
              vector<int>::const_iterator q,
              int ndim, 
              const int* dims) {
  for(int i=0; i<ndim; ++i) {
    int k=p[i]+q[i];
    if(k<0 || k>=dims[i]) {
      return false;
    }
  }
  return true;
}

vector<int>
makeNeighborhood10(int ndim, const int* dims) {
  if(ndim<3) {
    //mexPrintf("makeNeighborhood10: needs 3 dimensions.\n");
    vector<int> empty;
    return empty;
  }

  const int n=10;
  int id[n]={-dims[0]*dims[1],-dims[0]-1,-dims[0],-dims[0]+1,\
    -1,1,dims[0]-1,dims[0],dims[0]+1,dims[0]*dims[1]};

  vector<int> nbh;
  for(int i=0; i<n; ++i)
    nbh.push_back(id[i]);

  return nbh;
}

/*
Generalization of 4-neighborhood to N-dimensions
*/
vector<int>
MakeFourNeighborhood(int ndim,
                     const int* dims) {

  vector<int> vindex;
  for(int i=0; i<ndim; ++i) {
    vector<int> vsub(ndim,0);
    vsub[i]=-1;
    vindex.push_back(Sub2IndCentered(vsub,ndim,dims));
    vsub[i]=1;
    vindex.push_back(Sub2IndCentered(vsub,ndim,dims));
  }

  return vindex;
}

/*
Generalization of causal 4-neighborhood to N-dimensions
*/
vector<int>
MakeCausalFourNeighborhood(int ndim,
                     const int* dims) {

  vector<int> vindex;
  for(int i=0; i<ndim; ++i) {
    vector<int> vsub(ndim,0);
    vsub[i]=-1;
    vindex.push_back(Sub2IndCentered(vsub,ndim,dims));
  }

  return vindex;
}

vector<int>
MakeCausalNeighborhood(int ndim,
                       const int* dims) {
  int num_comb = (Round(pow(3.0, (double)ndim)-1.0))/2;
  vector<int> vindex;
  for(int n=0; n<num_comb; ++n) {
    vector<int> vsub(ndim,0);
    int m=n;
    for(int k=0; k<ndim; ++k) {
      int rem = m % 3;
      if(rem==0)
        vsub[k]=-1;
      else if(rem==1)
        vsub[k]=0;
      else 
        vsub[k]=1;
      m = m / 3;
	  //mexPrintf("MakeCausalNeighborhood: n = %d, dim = %d, rem = %d, sub = %d\n", n, k, rem, vsub[k]);
    }
    int index=Sub2IndCentered(vsub,ndim,dims);
    vindex.push_back(index);
  }

  for(int i=0; i<vindex.size(); ++i) {
    //mexPrintf("MakeCausalNeighborhood: %d - ", vindex[i]);
    vector<int> vsub=Ind2SubCentered(vindex[i],ndim,dims);
    //for(int j=0; j<vsub.size(); ++j)
    //  mexPrintf("%d ", vsub[j]);
    //mexPrintf("\n");
  }

  return vindex;
}

vector<int>
MakeAntiCausalNeighborhood(int ndim,
                           const int* dims) {
  int num_comb = (Round(pow(3.0, (double)ndim)-1.0))/2;
  vector<int> vindex;
  for(int n=0; n<num_comb; ++n) {
    vector<int> vsub(ndim,0);
    int m=n;
    for(int k=0; k<ndim; ++k) {
      int rem = m % 3;
      if(rem==0)
        vsub[k]=1;
      else if(rem==1)
        vsub[k]=0;
      else 
        vsub[k]=-1;
      m = m / 3;
    }
    int index=Sub2IndCentered(vsub,ndim,dims);
    vindex.push_back(index);
  }

  for(int i=0; i<vindex.size(); ++i) {
    //printf("MakeAntiCausalNeighborhood: %d\n", vindex[i]);
  }

  return vindex;
}

/*
Generalization of 8-neighborhood to N-dimensions
*/
vector<int>
MakeEightNeighborhood(int ndim,
                     const int* dims) {

  vector<int> vindex=MakeCausalNeighborhood(ndim,dims);
  vector<int> vindex2=MakeAntiCausalNeighborhood(ndim,dims);
  vindex.insert(vindex.end(),vindex2.begin(),vindex2.end());

  //for(int i=0; i<vindex.size(); ++i) {
  //  mexPrintf("MakeEightNeighborhood: %d\n", vindex[i]);
  //}

  return vindex;
}

/*
Generalization of 9-neighborhood to N-dimensions
*/
vector<int>
MakeNineNeighborhood(int ndim,
                     const int* dims) {

  vector<int> vindex=MakeCausalNeighborhood(ndim,dims);
  vector<int> vindex2=MakeAntiCausalNeighborhood(ndim,dims);
  vindex.insert(vindex.end(),0); 
  vindex.insert(vindex.end(),vindex2.begin(),vindex2.end());

  return vindex;
}

/*
rectanglular neighborhood of an arbitrary size
*/
vector<int>
MakeNeighborhood(const int* width,
				 int ndim,
				 const int* dims) 
{
	int i, j;
	vector<int> fwidth(ndim);
	int n = 1;
	for(i=0; i<ndim; ++i)
	{
		fwidth[i] = 2*width[i]+1;
		n *= fwidth[i];
	}
	vector<int> vindex;
	vector<int> vsub(ndim);
	for(i=0; i<n; ++i)
	{
		int m = i;
		bool bSelf = true;
		for(j=0; j<ndim; ++j)
		{
			vsub[j] = m % fwidth[j] - width[j];
			m /= fwidth[j];
			if(vsub[j] != 0)
			{
				bSelf = false;
			}
		}
		if(bSelf)
		{
			//self is not its neighbor
			continue;
		}

		int index=Sub2IndCentered(vsub,ndim,dims);
		if(index != 0)
		{
			vindex.push_back(index);
		}
	}
	return vindex;
}

vector<int>
MakeNeighborhood(int width,
				 int ndim,
				 const int* dims) 
{
	int fwidth = 2*width + 1;
	int n = (int)(pow((double)fwidth, (double)ndim));
	vector<int> vindex;
	vector<int> vsub(ndim);
	int i, j;
	for(i=0; i<n; ++i)
	{
		int m = i;
		bool bSelf = true;
		for(j=0; j<ndim; ++j)
		{
			vsub[j] = m % fwidth - width;
			m /= fwidth;
			if(vsub[j] != 0)
			{
				bSelf = false;
			}
		}
		if(bSelf)
		{
			//self is not its neighbor
			continue;
		}

		int index=Sub2IndCentered(vsub,ndim,dims);
		if(index != 0)
		{
			vindex.push_back(index);
		}
	}
	return vindex;
}

