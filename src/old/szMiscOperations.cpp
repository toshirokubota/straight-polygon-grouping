
#include <mex.h>
#include <math.h>
#include <algorithm>
using namespace std;

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>

int GetOffsetX6(int n) 
{
  switch(n) {
  case 1:
    return -1;
  case 2:
    return 1;
  case 0:
  case 3:
  case 4:
  case 5:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetY6(int n) 
{
  switch(n) {
  case 0:
    return -1;
  case 3:
    return 1;
  case 1:
  case 2:
  case 4:
  case 5:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetZ6(int n) 
{
  switch(n) {
  case 4:
    return -1;
  case 5:
    return 1;
  case 0:
  case 1:
  case 2:
  case 3:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
  //const int zoff[6]={0,0,0,0,-1,1};
  //return zoff[n];
}

int GetOffsetX10(int n) 
{
  switch(n) {
  case 1:
  case 4:
  case 6:
    return -1;
  case 3:
  case 5:
  case 8:
    return 1;
  case 0:
  case 2:
  case 7:
  case 9:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetY10(int n) 
{
  switch(n) {
  case 1:
  case 2:
  case 3:
    return -1;
  case 6:
  case 7:
  case 8:
    return 1;
  case 0:
  case 4:
  case 5:
  case 9:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetZ10(int n) 
{
  switch(n) {
  case 0:
    return -1;
  case 9:
    return 1;
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
  //const int zoff[6]={0,0,0,0,-1,1};
  //return zoff[n];
}

int GetOffsetX26(int n) 
{
  if(n<9)
    return (n%3)-1;
  else if(n<17) {
    switch(n) {
    case 9:
    case 12:
    case 14:
      return -1;
    case 10:
    case 15:
      return 0;
    case 11:
    case 13:
    case 16:
      return 1;
	default:
		return 0;
    }
  }
  else 
    return (n-17)%3 - 1;
}
 
int GetOffsetY26(int n) 
{
  if(n<9)
    return (n/3)-1;
  else if(n<17) {
    switch(n) {
    case 9:
    case 10:
    case 11:
      return -1;
    case 12:
    case 13:
      return 0;
    case 14:
    case 15:
    case 16:
      return 1;
	default:
		return 0;
    }
  }
  else 
    return (n-17)/3-1;
}
 
int GetOffsetZ26(int n) {
  if(n<9)
    return -1;
  else if(n<17)
    return 0;
  else //if(n<26)
    return 1;
}

vector<real> 
SymmetricExtension(const vector<real>& A, 
				   const int* dims)
{
	int dims2[]={dims[0]+2, dims[1]+2};
	vector<real> B(dims2[0]*dims2[1], 0);
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			SetData2(B, x+1, y+1, dims2[0], dims2[1], GetData2(A, x, y, dims[0], dims[1], (real)0));
		}
	}
	for(int y=1; y<dims2[1]-1; ++y)
	{
		real a = GetData2(B, 1, y, dims2[0], dims2[1], (real)0);
		real b = GetData2(B, 2, y, dims2[0], dims2[1], (real)0);
		SetData2(B, 0, y, dims2[0], dims2[1], 2*a-b);
		real c = GetData2(B, dims2[0]-2, y, dims2[0], dims2[1], (real)0);
		real d = GetData2(B, dims2[0]-3, y, dims2[0], dims2[1], (real)0);
		SetData2(B, dims2[0]-1, y, dims2[0], dims2[1], 2*c-d);
	}
	for(int x=0; x<dims2[0]; ++x)
	{
		real a = GetData2(B, x, 1, dims2[0], dims2[1], (real)0);
		real b = GetData2(B, x, 2, dims2[0], dims2[1], (real)0);
		SetData2(B, x, 0, dims2[0], dims2[1], 2*a-b);
		real c = GetData2(B, x, dims2[1]-2, dims2[0], dims2[1], (real)0);
		real d = GetData2(B, x, dims2[1]-3, dims2[0], dims2[1], (real)0);
		SetData2(B, x, dims2[1]-1, dims2[0], dims2[1], 2*c-d);
	}
	return B;
}

void
CurvatureImage(vector<real>& K,
			   const vector<real>& U,
			   const int* dims) 
{
	vector<real> V = SymmetricExtension(U, dims);
	int dims2[] = {dims[0]+2, dims[1]+2};
	for(int y1=0; y1<dims[1]; ++y1)
	{
		for(int x1=0; x1<dims[0]; ++x1)
		{
			int x = x1 + 1;
			int y = y1 + 1;
			real c = GetData2(V, x, y, dims2[0], dims2[1], (real)0);
			real dx = (GetData2(V, x+1, y, dims2[0], dims2[1], c) - GetData2(V, x-1, y, dims2[0], dims2[1], c))/2.0; ///2.0;
			real dy = (GetData2(V, x, y+1, dims2[0], dims2[1], c) - GetData2(V, x, y-1, dims2[0], dims2[1], c))/2.0; ///2.0;
			real dxx = (GetData2(V, x+1, y, dims2[0], dims2[1], c) + GetData2(V, x-1, y, dims2[0], dims2[1], c) - 2*c); ///4.0;
			real dyy = (GetData2(V, x, y+1, dims2[0], dims2[1], c) + GetData2(V, x, y-1, dims2[0], dims2[1], c) - 2*c); ///4.0;
			real dxy = (GetData2(V, x-1, y-1, dims2[0], dims2[1], c) + GetData2(V, x+1, y+1, dims2[0], dims2[1], c)\
				- GetData2(V, x+1, y-1, dims2[0], dims2[1], c) - GetData2(V, x-1, y+1, dims2[0], dims2[1], c))/4.0;
			real hyp=sqrt(dx*dx+dy*dy); //add one to the denominator
			real cv=dyy*(dx*dx) + dxx*(dy*dy) - 2*dx*dy*dxy; //add one to dx^2 and dy^2 terms
			if(hyp < 1.0e-6)
			{
				cv = 0;
			}
			else
			{
				cv/=(hyp*hyp*hyp);
			}
			SetData2(K, x1, y1, dims[0], dims[1], cv);
		}
	}
}

/*float
Distance(const CParticleF& p1, const CParticleF& p2)
{
	return sqrt((float)(p1.m_X-p2.m_X)*(p1.m_X-p2.m_X)+(p1.m_Y-p2.m_Y)*(p1.m_Y-p2.m_Y));
}*/

void
RemoveJunctions(vector<unsigned char>& I,
				const int* dims)
{
	int xoff[]={-1, 0, 1, 1, 1, 0, -1, -1};
	int yoff[]={-1, -1, -1, 0, 1, 1, 1, 0};
	int numChanged;
	do
	{
		vector<unsigned char> E = I;
		numChanged = 0;
		for(int y=0; y<dims[1]; ++y)
		{
			for(int x=0; x<dims[0]; ++x)
			{
				if(GetData2(E, x, y, dims[0], dims[1], (unsigned char)0) == 0)
				{
					continue;
				}
				unsigned char prev = GetData2(E, x+xoff[7], y+yoff[7], dims[0], dims[1], (unsigned char)0);
				int cnt = 0;
				for(int k=0; k<8; ++k)
				{
					unsigned char cval = GetData2(E, x+xoff[k], y+yoff[k], dims[0], dims[1], (unsigned char)0);
					if(prev==0 && cval)
					{
						cnt++;
					}
					prev = cval;
				}
				if(cnt > 2)
				{
					SetData2(I, x, y, dims[0], dims[1], (unsigned char)0);
					numChanged++;
				}
			}
		}
	}
	while(numChanged>0);
}

