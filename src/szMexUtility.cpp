#ifdef MEX_DLL
#include <mex.h>
#endif
#include <szMexUtility.h>

int 
Mod(int a, int b) {
  if(a>0)
    return a % b;
  else {
    int c=(-a) % b;
    if(c==0)
      return 0;
    else
      return b-c;
  }
}

int
numberOfElements(int ndim, const int* dims) {
  if(ndim<=0 || dims==NULL)
    return 0;

  int n=1;
  for(int i=0; i<ndim; ++i)
    n*=dims[i];

  return n;
}

int
Sub2Ind(int x, int y, int z, const int* dims) {
  return z*dims[0]*dims[1]+y*dims[0]+x;
}

void
Ind2Sub(int& x, int& y, int& z, int ind, const int* dims) {

  x=ind % dims[0];
  ind/=dims[0];
  y=ind % dims[1];
  z=ind/dims[1];

  return;
}

int
Sub2Ind(const vector<int>& vsub, int ndim, const int* dims) {
  int ind=0;
  int stride=1;
  for(int i=0; i<vsub.size(); ++i) {
    ind+=vsub[i]*stride;
    stride*=dims[i];
  }
  return ind;
}

vector<int>
Ind2Sub(int ind, int ndim, const int* dims) {

  vector<int> vsub(ndim,0);
  for(int i=0; i<ndim; ++i) {
    vsub[i] = ind % dims[i];
    ind /= dims[i];
  }
  if(ind)
    vsub.push_back(ind); //overflow

  return vsub;
}

void
Ind2Sub(vector<int>& vsub, int ind, int ndim, const int* dims) 
{

	for(int i=0; i<ndim; ++i) 
	{
		if(vsub.size()< i)
			vsub.push_back(ind % dims[i]);
		else
			vsub[i] = ind % dims[i];
		ind /= dims[i];
	}
	if(ind)
		vsub.push_back(ind); //overflow
}

/*
This sub2ind version computes subscript with respect to the center of the N-d data.
It treats the mid point of the index as the origin ([0,0,...,0]).
*/
int
Sub2IndCentered(const vector<int>& vsub, int ndim, const int* dims) {
  vector<int> vsub2(ndim);
  vector<int> vsub3(ndim);
  for(int i=0; i<ndim; ++i) {
    vsub2[i] = vsub[i] + dims[i]/2;
    vsub3[i] = dims[i]/2;
  }
  int ind = Sub2Ind(vsub2, ndim, dims);
  int ind2 = Sub2Ind(vsub3, ndim, dims);
  ind-=ind2;

  /*mexPrintf("Sub2IndCentered: (");
  for(i=0; i<ndim; ++i) 
    mexPrintf("%d ", vsub[i]);
  mexPrintf(") -> %d\n", ind);*/

  return ind;
}

/*
This ind2sub version computes subscript with respect to the center of the N-d data.
It treats the mid point of the index as the origin ([0,0,...,0]).
*/
vector<int>
Ind2SubCentered(int ind, int ndim, const int* dims) {
  vector<int> vsub2(ndim);
  int i;
  for(i=0; i<ndim; ++i) {
    vsub2[i] = dims[i]/2;
  }
  int off=Sub2Ind(vsub2,ndim,dims);

  ind+=off;
  vector<int> vsub = Ind2Sub(ind, ndim, dims);
  for(i=0; i<ndim; ++i) {
    vsub[i] = vsub[i] - dims[i]/2;
  }

  return vsub;
}

/*
compute the dimensions of the volume when voxel size is to be changed.
it is used during volume interpolation process.
*/
void
getNewDimension(int* newdims,
				vector<float>& vsize_new,
				vector<float>& vsize_old,
				int ndim,
				const int* olddims)
{
	int i;
	for(i=0; i<ndim; ++i)
	{
		if(olddims[i] % 2)
		{
			//keep the center location unchanged
			float r = ((olddims[i]-1)/2) * vsize_old[i] / vsize_new[i];
			newdims[i] = (int)r * 2 + 1;
			//printf("%d getNewDimension: odd\n", i);
		}
		else
		{  //keep the location one location off unchanged
			float r1 = (olddims[i]/2 - 1) * vsize_old[i] / vsize_new[i];
			float r2 = (olddims[i]/2 + 1) * vsize_old[i] / vsize_new[i];
			newdims[i] = (int)r1 + (int)r2 + 1;
			//printf("%d getNewDimension: even\n", i);
		}
	}
}

double
rndm(int jd)
{
   /*static int mdig = 16;*/	/* # digits used for integers */
/*
   M1= 2**(MDIG-2) + (2**(MDIG-2)-1),  M2 = 2**(MDIG/2)
*/

   static int m[17] =
   {30788, 23052, 2053, 19346, 10646, 19427, 23975,
   19049, 10949, 19693, 29746, 26748, 2796, 23890, 29168, 31924, 16499};

   static int i = 4, j = 16, m1 = 32767, m2 = 256;
   static double invm1 = .00003051850947599719;  /* 1.0 / m1 */

   int k, jseed, k0, k1, j0, j1;

   if (jd != 0)
   {
      jseed = (jd < 0) ? -jd : jd;
      jseed = (jseed < m1) ? jseed : m1;
      if (jseed % 2 == 0)
	 --jseed;
      k0 = 9069 % m2;
      k1 = 9069 / m2;
      j0 = jseed % m2;
      j1 = jseed / m2;
#if 0
      printf("jseed= %d, k0= %d, k1= %d, j0= %d, j1= %d\n", jseed, k0, k1, j0, j1);
#endif
      for (i = 0; i < 17; ++i)
      {
	 jseed = j0 * k0;
	 j1 = ((jseed / m2) + j0 * k1 + j1 * k0) % (m2 / 2);
	 j0 = jseed % m2;
	 m[i] = j0 + m2 * j1;
#if 0
	 printf("i= %d, m[i]= %d\n", i, m[i]);
#endif
      }
      i = 4;
      j = 16;
   }
   k = m[i] - m[j];
   if (k < 0)
      k += m1;
   m[j] = k;
   if (--i < 0)
      i = 16;
   if (--j < 0)
      j = 16;
   return (double) k *  invm1;
}
