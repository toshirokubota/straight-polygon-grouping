#include <NonMaximumSuppression.h>
using namespace std;

vector<float>
CannyNonMaximumSuppression(const vector<float>& dmag,
                           const vector<float>& dy, 
						   const vector<float>& dx,
						   const int* dims) 
{
	vector<float> res(dims[0]*dims[1],0);
	int i,j;
	for (i=1; i<dims[1]-1;++i)
	{
		for (j=1; j<dims[0]-1; ++j)
		{
			float ux=GetData2(dx,i, j, dims[0], dims[1], 0.0f);
			float uy=GetData2(dy,i, j, dims[0], dims[1], 0.0f);
			float or=atan2(uy,ux);

			float gc=GetData2(dmag,i, j, dims[0], dims[1], 0.0f);
			float gn=GetData2(dmag,i-1, j, dims[0], dims[1], 0.0f);
			float gs=GetData2(dmag,i+1, j, dims[0], dims[1], 0.0f);
			float gw=GetData2(dmag,i, j-1, dims[0], dims[1], 0.0f);
			float ge=GetData2(dmag,i, j+1, dims[0], dims[1], 0.0f);
			float gne=GetData2(dmag,i-1, j+1, dims[0], dims[1], 0.0f);
			float gse=GetData2(dmag,i+1, j+1, dims[0], dims[1], 0.0f);
			float gnw=GetData2(dmag,i-1, j-1, dims[0], dims[1], 0.0f);
			float gsw=GetData2(dmag,i+1, j-1, dims[0], dims[1], 0.0f);
			float g0;
			if(ux*uy>0) 
			{
				if(Abs(ux)<Abs(uy)) 
				{
					if((g0=Abs(uy*gc))<Abs(ux*gse+(uy-ux)*gs) ||
						g0<Abs(ux*gnw+(uy-ux)*gn)) //used to be <=
						continue;
				}
				else 
				{
					if((g0=Abs(ux*gc))<Abs(uy*gse+(ux-uy)*ge) ||
						g0<Abs(uy*gnw+(ux-uy)*gw)) //used to be <=
						continue;
				}
			}
			else 
			{
				if(Abs(ux)<Abs(uy)) 
				{
					if((g0=Abs(uy*gc))<Abs(ux*gne-(uy+ux)*gn) ||
						g0<Abs(ux*gsw-(uy+ux)*gs)) //used to be <=
						continue;
				}
				else 
				{
					if((g0=Abs(ux*gc))<Abs(uy*gne-(ux+uy)*ge) ||
						g0<Abs(uy*gsw-(ux+uy)*gw)) //used to be <=
						continue;
				}
			}
			SetData2(res, i, j, dims[0], dims[1], gc); 
		}
	} 
	return res;
}

float
ThresholdSelection(const vector<float> grd, 
				   float percent,
				   const int* dims) 
{
	int i;
	int nump=dims[0] * dims[1];
	float maxval = grd[0];
	float minval= maxval;
	for(i=1; i<nump; ++i) 
	{
		float val=grd[i];
		if(maxval<val)
			maxval=val;
		if(minval>val)
			minval=val;
	}
	if(maxval==minval) //constant image
		return maxval;

	int numbin=256;
	vector<float> vhist(numbin,0);
	for(i=0; i<nump; ++i)
	{
		float val=grd[i];
		int id=(int)((numbin-1)*(val-minval)/(maxval-minval));
		vhist[id]++;
	}
	bool found=false;
	int idx=numbin;
	do 
	{
		int target=Round(nump*percent);
		int sum=0;
		for(i=0; i<numbin; ++i) 
		{
			sum+=vhist[i];
			if(sum>target)
				break;
		}
		//for image with small number of edges, (in particular synthetic one)
		//the selection may still be the same with minval.
		//if so, increase the percentile and try again.
		if(i>0) 
		{
			found=true;
			idx = i;
		}
		else
			percent=1-(1-percent)/2;
	}
	while(!found);
	return (float)idx*(maxval-minval)/numbin+minval;
}

void
follow(vector<float>& trace, 
	   const vector<float>& edge, 
	   vector<unsigned char>& flag, 
	   int y, int x, 
	   float low,
	   const int* dims) 
{
	SetData2(flag, x, y, dims[0], dims[1], (unsigned char)1);
	SetData2(trace, x, y, dims[0], dims[1], GetData2(edge, x, y, dims[0], dims[1], 0.0f));
	for(int i=y-1; i<=y+1; ++i) 
	{
		for(int j=x-1; j<=x+1; ++j) 
		{
			if(!GetData2(flag, j, i, dims[0], dims[1], (unsigned char)1))
			{
				if(GetData2(edge, j, i, dims[0], dims[1], (float)0)>low)
				{
					follow(trace, edge, flag, i, j, low, dims);
				}
			}
		}
	}
}

vector<float>
EdgeTrace(const vector<float>& edge, 
		  float high, float low,
		  const int* dims) 
{
	vector<float> res(dims[0]*dims[1],.0f);
	vector<unsigned char> flag(dims[0]*dims[1],(unsigned char)0);
	for(int i=0; i<dims[1]; ++i) 
	{
		for(int j=0; j<dims[0]; ++j) 
		{
			if(!GetData2(flag, j, i, dims[0], dims[1], (unsigned char)1) && GetData2(edge, j, i, dims[0], dims[1], (float)0) > high)
			{
				follow(res, edge, flag, i, j, low, dims);
			}
		}
	}
	return res;
}


/*
Morphological thinning based on 
Gua & Hall, CACM, 1989, p359-373.
One operation is enough after the non-maximum suppression.
*/ 
vector<unsigned char>
CannyThinning(const vector<unsigned char>& image,
			  const int* dims) 
{
	const int NumNeighbors=8;
	int ax[NumNeighbors]={1,1,0,-1,-1,-1,0,1};
	int ay[NumNeighbors]={0,-1,-1,-1,0,1,1,1};

	vector<unsigned char> res=image;
	int i,j,m;
	for(int pass=0; pass<2; ++pass) 
	{
		for(i=1; i<dims[1]-1; ++i) 
		{
			for(j=1; j<dims[0]-1; ++j) 
			{
				if(GetData2(res, j, i, dims[0], dims[1], (unsigned char)0))
				{
					int cp=0;
					int np1=0,np2=0;
					for(m=1; m<=NumNeighbors/2; m++) {
						int y1=i+ay[(2*m-2)%NumNeighbors];
						int x1=j+ax[(2*m-2)%NumNeighbors];
						int y2=i+ay[(2*m-1)%NumNeighbors];
						int x2=j+ax[(2*m-1)%NumNeighbors];
						int y3=i+ay[(2*m)%NumNeighbors];
						int x3=j+ax[(2*m)%NumNeighbors];
						unsigned char p1=GetData2(res, x1, y1, dims[0], dims[1], (unsigned char)0);
						unsigned char p2=GetData2(res, x2, y2, dims[0], dims[1], (unsigned char)0);
						unsigned char p3=GetData2(res, x3, y3, dims[0], dims[1], (unsigned char)0);
						if(!p1 && (p2||p3))
							cp++;
						if(p1||p2)
							np1++;
						if(p2||p3)
							np2++;
					}
					int np=Min(np1,np2);
					if(cp==1 && np>=2 && np<=3) 
					{
						if(pass==0) 
						{
							unsigned char p2=GetData2(res, j+ax[1], i+ay[1], dims[0], dims[1], (unsigned char)0);
							unsigned char p3=GetData2(res, j+ax[2], i+ay[2], dims[0], dims[1], (unsigned char)0);
							unsigned char p8=GetData2(res, j+ax[7], i+ay[7], dims[0], dims[1], (unsigned char)0);
							unsigned char p1=GetData2(res, j+ax[0], i+ay[0], dims[0], dims[1], (unsigned char)0);
							if(!((p2 || p3 || !p8) && p1))
							{
								SetData2(res, j, i, dims[0], dims[1], (unsigned char)0);
							}
						}
						else if(pass==1) 
						{
							unsigned char p4=GetData2(res, j+ax[3], i+ay[3], dims[0], dims[1], (unsigned char)0);
							unsigned char p5=GetData2(res, j+ax[4], i+ay[4], dims[0], dims[1], (unsigned char)0);
							unsigned char p6=GetData2(res, j+ax[5], i+ay[5], dims[0], dims[1], (unsigned char)0);
							unsigned char p7=GetData2(res, j+ax[6], i+ay[6], dims[0], dims[1], (unsigned char)0);
							if(!((p6 || p7 || !p4) && p5))
							{
								SetData2(res, j, i, dims[0], dims[1], (unsigned char)0);
							}
						}                 
					}
				}
			}
		}
	}
	return res;
}

#include <mex.h>

vector<unsigned char>
thinning(const vector<float>& im,
		 const vector<float>& dx,
		 const vector<float>& dy,
		 const int* dims)
{
	float percent = 0.7f;
	float ratio = 0.4f;
	vector<float> emag=CannyNonMaximumSuppression(im, dy, dx, dims);
	float high=ThresholdSelection(im, percent, dims);
	float low=high*ratio;
	printf("low thres = %f, high thres = %f\n", low, high);
	emag=EdgeTrace(emag, high, low, dims);

	vector<unsigned char> edge(dims[0]*dims[1], (unsigned char)0);
	for(int i=0; i<dims[0]*dims[1];++i)
	{
		edge[i] = emag[i] > 0 ? 1: 0;
	}

	edge=CannyThinning(edge, dims);
	edge=CannyThinning(edge, dims);
	return edge;
}
