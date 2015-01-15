#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
using namespace std;
#include <szParticle.h>
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConnectedComponentC.h>
#include <szIndexedData.h>
#include <szContour.h>

const int XOffset8R[] = {-1, 0, 1, 1, 1, 0, -1, -1};
const int YOffset8R[] = {-1, -1, -1, 0, 1, 1, 1, 0};

//const int  NumNeighbors8 = 8;
enum PointType {None, Normal, End, Junction, Start, Finish};

/*
Priopritize directions based on the current position in relation to the previous one.
It is used to guide tracing a contour.
*/
struct TraceIterator
{
	TraceIterator(int x, int y, int px, int py)
	{
		int dx = x - px;
		int dy = y - py;
		int ix = dx > 0 ? 1: (dx < 0 ? -1: 0);
		int iy = dy > 0 ? 1: (dy < 0 ? -1: 0);
		if(ix == 0 && iy==0) index = 0;
		else
		{
			if(ix==-1 && iy==-1) index = 0;
			else if(ix==0 && iy==-1) index = 1;
			else if(ix==1 && iy==-1) index = 2;
			else if(ix==-1 && iy==0) index = 7;
			else if(ix==1 && iy==0) index = 3;
			else if(ix==-1 && iy==1) index = 6;
			else if(ix==0 && iy==1) index = 5;
			else if(ix==1 && iy==1) index = 4;
		}
		count = 0;
	}
	TraceIterator(int x, int y)
	{
		*this = TraceIterator(x, y, x, y);
	}
	TraceIterator(const CParticleF& q, const CParticleF& p)
	{
		*this = TraceIterator(q.m_X, q.m_Y, p.m_X, p.m_Y);
	}
	bool Next(int& x, int& y)
	{
		if(count >= 8) return false; 
		else
		{
			int k;
			if(count % 2 == 1)
			{
				k = (index - (count+1)/2 + 8) % 8;
			}
			else
			{
				k = (index + count/2) % 8;
			}
			x = XOffset8R[k];
			y = YOffset8R[k];
			count ++;
			return true;
		}
	}
private:
	int index;
	int count;
};

/*
For each edge point, go around its 8-neighbors and see if there are three pixels in series that 
are edges. Then, this poit can be removed.
*/ 
void
Thinning(vector<PointType>& J, const int* dims) 
{
	for(int i=0; i<dims[1]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			if(GetData2(J, j, i, dims[0], dims[1], None) != None)
			{
				PointType vals[8];
				for(int k=0; k< NumNeighbors8; ++k)
				{
					int i2 = i+YOffset8R[k];
					int j2 = j+XOffset8R[k];
					vals[k] = GetData2(J, j2, i2, dims[0], dims[1], None);
				}
				int first = -1;
				for(int k=0; k< NumNeighbors8; ++k)
				{
					if(vals[k] != None)
					{
						first = k;
						break;
					}
				}
				if(first < 0) continue;
				if(vals[1] && vals[3] || vals[3] && vals[5] || vals[5] && vals[7] || vals[7] && vals[1])
				{
					SetData2(J, j, i, dims[0], dims[1], None);
					continue;
				}
				int maxseq=0;
				int count = 0;
				for(int k=0; k< NumNeighbors8; ++k)
				{
					int k2 = (first + k) % NumNeighbors8;
					if(vals[k2] != None)
					{
						count++;
						if(count > maxseq) maxseq = count;
					}
					else
					{
						count=0;
					}
				}

				if(maxseq >= 3) //this point can be removed.
				{
					SetData2(J, j, i, dims[0], dims[1], None);
				}
			}
		}
	}
}

/*
Set Junction, End, and Normal PointType for each edge pixel.
*/
vector<PointType>
SetJunctionEndFlags(const vector<PointType>& I,
					const int* dims)
{
	vector<PointType> J(dims[0] * dims[1], None);
	for(int i=0; i<dims[1]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			if(GetData2(I, j, i, dims[0], dims[1], None) != None)
			{
				int count = 0;
				int count2 = 0;
				PointType val1 = GetData2(I, j+XOffset[0], i+YOffset[0], dims[0], dims[1], None);
				//PointType vals[8]={None,None,None,None,None,None,None,None};
				for(int k=1; k<= NumNeighbors8; ++k)
				{
					int i2 = i+YOffset8R[k % NumNeighbors8];
					int j2 = j+XOffset8R[k % NumNeighbors8];
					PointType val2 = GetData2(I, j2, i2, dims[0], dims[1], None);
					if(val1==None && val2!=None)
					{
						count++;
					}
					if(val2!=None)
					{
						count2++;
					}
					val1 = val2;
				}
				if(count > 2)
				{
					SetData2(J, j, i, dims[0], dims[1], Junction);
				}
				else if(count == 1)
				{
					SetData2(J, j, i, dims[0], dims[1], End);
				}
				else if(count2 == 0)
				{
					SetData2(J, j, i, dims[0], dims[1], None);
				}
				else
				{
					SetData2(J, j, i, dims[0], dims[1], Normal);
				}
			}
		}
	}
	return J;
}

void
trace(vector<PointType>& J,
	  Contour& c,
	  const CParticleF& pnt,
	  const int* dims)
{
	CParticleF prev = c.points.size()>0 ? c.points[c.points.size()-1]: pnt;
	c.points.push_back(pnt);
	PointType t = GetData2(J, pnt.m_X, pnt.m_Y, dims[0], dims[1], Normal);

	if(t==End || t==Junction)
	{
		return;
	}
	else if(t==Start)
	{
		SetData2(J, pnt.m_X, pnt.m_Y, dims[0], dims[1], Finish);
	}

	TraceIterator it(pnt, prev);
	int offx, offy;
	while(it.Next(offx, offy))
	{
		int i2 = pnt.m_Y+offy;
		int j2 = pnt.m_X+offx;
		if(GetData2(J, j2, i2, dims[0], dims[1], None)!=None)
		{
			CParticleF q(j2, i2);
			vector<CParticleF>::iterator it = find(c.points.begin(), c.points.end(), q);
			if(it==c.points.end())
				///TK! Do we need this condition?-> || (GetData2(J, j2, i2, dims[0], dims[1], Normal)==Finish && c.points.size()>3))
			{
				trace(J, c, q, dims);
				break;
			}
		}
	}
	return;
}

Contour
traceContour(vector<PointType>& J,
			 const CParticleF& p,
			 const int* dims)
{
	Contour c;
	PointType t = GetData2(J, p.m_X, p.m_Y, dims[0], dims[1], Normal);
	SetData2(J, p.m_X, p.m_Y, dims[0], dims[1], Start);
	trace(J, c, p, dims);
	SetData2(J, p.m_X, p.m_Y, dims[0], dims[1], t);

	return c;
}

template <class Pr>
void
ProcessIt(vector<PointType>& J,
		  vector<Contour>& contours,
		  int minlen,
		  Pr pred,       //Use () operator and select initial points (in terms of the PointType) for trace.
		  const int* dims)
{
	int count = 0;
	//trace from end points.
	while(true)
	{
		//SetJunctionEndFlags(J, dims);
		//if(count==0) J0 = J; //save the original
		bool bChanged = false;
		for(int i=0; i<dims[1]; ++i)
		{
			for(int j=0; j<dims[0]; ++j)
			{
				if(pred(GetData2(J, j, i, dims[0], dims[1], Normal)))
				{
					CParticleF p(j, i, 0);
					Contour c = traceContour(J, p, dims);
					if(c.points.size()>=minlen)
					{
						contours.push_back(c);
						bChanged = true;
					}
					for(int m=0; m<c.points.size(); ++m)
					{
						int x = c.points[m].m_X;
						int y = c.points[m].m_Y;
						PointType type = GetData2(J, x, y, dims[0], dims[1], None);
						if(type != Junction)
						{
							SetData2(J, x, y, dims[0], dims[1], None);
						}
					}
				}
			}
		}
		count++;
		if(bChanged == false) break;
	}
}

struct EndPredicator
{
	bool operator()(PointType type) { return type == End; }
};

struct EndJunctionPredicator
{
	bool operator()(PointType type) { return type == End || type == Junction; }
};

struct JunctionPredicator
{
	bool operator()(PointType type) {return type==Junction;}
};

struct OthersPredicator
{
	bool operator()(PointType type) {return type!=None;}
};

vector<CParticleF>
breakingPoints(Contour& contour, float thres)
{
	vector<CParticleF> pnts;
	if (contour.points.size() < 2)
	{
		return pnts;
	}
	CParticleF p1 = contour.points[0];
	CParticleF p2 = contour.points[contour.points.size() - 1];
	float maxd = 0.0f;
	int idx = 0;
	for (int i = 1; i < contour.points.size() - 1; i++)
	{
		float d = Distance2LineSegment(p1, p2, contour.points[i]);
		if (maxd < d)
		{
			idx = i;
			maxd = d;
		}
	}
	if (maxd > thres)
	{
		Contour a;
		for (int i = 0; i <= idx; ++i)
		{
			a.points.push_back(contour.points[i]);
		}
		vector<CParticleF> pa = breakingPoints(a, thres);
		Contour b;
		for (int i = idx; i < contour.points.size(); ++i)
		{
			b.points.push_back(contour.points[i]);
		}
		vector<CParticleF> pb = breakingPoints(b, thres);
		for (int i = 0; i < pa.size(); ++i)
		{
			pnts.push_back(pa[i]);
		}
		for (int i = 1; i < pb.size(); ++i)
		{
			pnts.push_back(pb[i]);
		}
	}
	else
	{
		pnts.push_back(p1);
		pnts.push_back(p2);
	}
	return pnts;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: C = ContourPrep(I)");
		return;
	}

	//Edge Image
	int ndimI;
	const int* dims;
	mxClassID classI;
	vector<unsigned char> I;
	LoadData(I,prhs[0],classI,ndimI,&dims);
	int minLen = 5;
	if(nrhs>=2) 
	{
		mxClassID classMode;
		ReadScalar(minLen,prhs[1],classMode);
	}
	float breakThres = 5.0;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(breakThres, prhs[2], classMode);
	}

	vector<Contour> contours;

	vector<PointType> J(I.size());
	for(int i=0; i<I.size(); ++i) J[i] = I[i] ? Normal: None;
	Thinning(J, dims);

	//first process from end points
	EndJunctionPredicator endPred;
	vector<PointType> K = SetJunctionEndFlags(J, dims);
	ProcessIt(K, contours, minLen, endPred, dims);

	//all we have left should be circular patterns.
	//first, process them from junctions
	//JunctionPredicator junctionPred;
	//vector<PointType> K2 = SetJunctionEndFlags(K, dims);
	//ProcessIt(K2, contours, minLen, junctionPred, dims);

	//all we have left should be circular patterns with no junction.
	//we can start trace from any point from such pattern.
	OthersPredicator othersPred;
	vector<PointType> K2 = SetJunctionEndFlags(K, dims);
	ProcessIt(K2, contours, minLen, othersPred, dims);

	vector<CParticleF> P;
	vector<int> E;
	for (int i = 0; i < contours.size(); ++i)
	{
		//if (i + 1 != 110) continue; ///TK. for debugging purpose
		vector<CParticleF> pnts = breakingPoints(contours[i], breakThres);
		vector<int> idx(pnts.size());
		for (int j = 0; j < pnts.size(); ++j)
		{
			vector<CParticleF>::iterator pi = find(P.begin(), P.end(), pnts[j]);
			if (pi == P.end())
			{
				P.push_back(pnts[j]);
				idx[j] = P.size() - 1;
			}
			else
			{
				idx[j] = distance(P.begin(), pi);
			}
		}
		for (int j = 1; j < pnts.size(); ++j)
		{
			E.push_back(idx[j - 1]);
			E.push_back(idx[j]);
		}
	}

	if (nlhs >= 1)
	{
		const int dims[] = { P.size(), 2 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], P[i].m_X);
			SetData2(F, i, 1, dims[0], dims[1], P[i].m_Y);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { E.size()/2, 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i<dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], E[2*i]+1);
			SetData2(F, i, 1, dims[0], dims[1], E[2*i+1]+1);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreContours(contours);
	}
	if (nlhs >= 4)
	{
		plhs[3]= StoreData(J, mxUINT8_CLASS, 2, dims);
	}
	if (nlhs >= 5)
	{
		plhs[4] = StoreData(K, mxUINT8_CLASS, 2, dims);
	}
	mexUnlock();
}

