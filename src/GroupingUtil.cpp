#include <GroupingUtil.h>
#include <algorithm>
using namespace std;
#include <szMexUtilityTemplate.h>
#include <szMiscOperations.h>
#include <LeastSquaresFitting.h>

real
Distance(real x, real y, real x2, real y2)
{
	return sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2));
}

int
findLeftActive(int k, 
			   const vector<EdgelQueue>& elist)
{
	Edgel* pe = elist[k].Top();
	for(int i=1; i<elist.size(); ++i)
	{
		int j = (k - i + elist.size()) % elist.size();
		Edgel* qe = elist[j].Top();
		if(qe && pe != qe && qe->Active)
		{
			return j;
		}
	}
	return k; //couldn't find.  Return the same index.
}

int
findRightActive(int k, 
			   const vector<EdgelQueue>& elist)
{
	Edgel* pe = elist[k].Top();
	for(int i=1; i<elist.size(); ++i)
	{
		int j = (k + i) % elist.size();
		Edgel* qe = elist[j].Top();
		if(qe && pe != qe && qe->Active)
		{
			return j;
		}
	}
	return k; //couldn't find.  Return the same index.
}

vector<CParticleF>
OrderEdgePoints(const vector<unsigned char>& I,
				const CParticleF& click,
				const int* dims)
{
	vector<CParticleF> vpoints;
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(GetData2(I, x, y, dims[0], dims[1], (unsigned char)0))
			{
				CParticleF p(x, y, 0);
				p.m_Life = Distance(p, click);
				vpoints.push_back(p);
			}
		}
	}
	sort(vpoints.begin(), vpoints.end());

	return vpoints;
}

vector<EdgelQueue>
BuildEdgeList(const vector<CParticleF>& vpoints,
			  const CParticleF& click,
			  const int* dims)
{
	vector<EdgelQueue> elist(AngleBinSize);
	EdgelFactory* factory = EdgelFactory::GetInstance();
	for(int i=0; i<vpoints.size(); ++i)
	{
		int x=vpoints[i].m_X;
		int y=vpoints[i].m_Y;
		if(x+1==35 && y+1==87)
			x+=0;
		real theta[4];
		theta[0] = atan2((real)(y-.5-click.m_Y), (real)(x-.5-click.m_X));
		theta[1] = atan2((real)(y-.5-click.m_Y), (real)(x+.5-click.m_X));
		theta[2] = atan2((real)(y+.5-click.m_Y), (real)(x-.5-click.m_X));
		theta[3] = atan2((real)(y+.5-click.m_Y), (real)(x+.5-click.m_X));
		int maxBin = (int)((theta[0] + PI)/(2*PI) * AngleBinSize);
		int minBin = maxBin;		
		for(int k=1; k<4; ++k)
		{
			int bin = (int)((theta[k] + PI)/(2*PI) * AngleBinSize);
			if(bin > maxBin) maxBin = bin;
			if(bin < minBin) minBin = bin;
		}
		Edgel* pe = factory->Make(x, y, 0, 1, 1, vpoints[i].m_Life);
		if(maxBin - minBin < AngleBinSize+minBin - maxBin)
		{
			for(int k=minBin; k<=maxBin; ++k)
			{
				elist[k].Enqueue(pe);
			}

		}
		else //go past the -pi and pi discontinuity 
		{
			int length = AngleBinSize+minBin - maxBin;
			for(int k=0; k<=length; ++k)
			{
				elist[(k + maxBin) % AngleBinSize].Enqueue(pe);
			}
		}
	}

	return elist;
}

/*
This version limits each edge to one bin.
*/
vector<EdgelQueue>
BuildEdgeListUnique(const vector<CParticleF>& vpoints,
					const CParticleF& click,
					const int* dims)
{
	vector<EdgelQueue> elist(AngleBinSize);
	EdgelFactory* factory = EdgelFactory::GetInstance();
	for(int i=0; i<vpoints.size(); ++i)
	{
		int x=vpoints[i].m_X;
		int y=vpoints[i].m_Y;
		real theta = atan2((real)(y-click.m_Y), (real)(x-click.m_X));
		int bin = (int)((theta + PI)/(2*PI) * AngleBinSize);
		bin = Min(bin, AngleBinSize-1);
		Edgel* pe = factory->Make(x, y, 0, 1, 1, vpoints[i].m_Life);
		elist[bin].Enqueue(pe);
	}

	return elist;
}

/*
Check if click is inside a polygon specified by edges.
The algorithm taken from http://alienryderflex.com/polygon/
*/
bool
enclosingClick(const vector<Edgel*> edges,
			   const CParticleF& click)
{
	bool oddNodes=false;
	real x = click.m_X;
	real y = click.m_Y;

	int j=edges.size()-1;
	for (int i=0; i<edges.size(); i++) 
	{
		real xi = edges[i]->X;
		real yi = edges[i]->Y;
		real xj = edges[j]->X;
		real yj = edges[j]->Y;
		if (yi<y && yj>=y || yj<y && yi>=y)
		{
			if (xi + (y-yi)/(yj-yi) * (xj-xi) < x)
			{
				oddNodes=!oddNodes; 
			}
		}
		j=i; 
	}

	return oddNodes; 
} 

real
computeElastica(const vector<Edgel*>& edges,
				bool& bFlag)
{
	if(edges.size()<=1)
	{
		bFlag = false;
		return 0;
	}
	real sum = 0;
	real alpha=10.0;
	for(int i=0; i<edges.size(); ++i)
	{
		int j = (i-1+edges.size()) % edges.size();
		int k = (i+1) % edges.size();
		Edgel* pe = edges[i];
		Edgel* left = edges[j];
		Edgel* right = edges[k];

		//compute the gaps
		real dL = Distance(pe->X, pe->Y, left->X, left->Y);
		real dR = Distance(pe->X, pe->Y, right->X, right->Y);

		//compute approximate curvature
		real cx = (left->X-pe->X)/(1+dL)+(right->X-pe->X)/(1+dR);
		real cy = (left->Y-pe->Y)/(1+dL)+(right->Y-pe->Y)/(1+dR);
		real cv = cx*cx + cy*cy;

		sum += dL + dR + alpha*cv;
	}
	bFlag = true;
	return sum;
}


real
computeCircularFitError(const vector<Edgel*>& edges,
						bool& bFlag)
{
	vector<real> A(edges.size()*3);
	vector<real> b(edges.size());
	for(int i=0; i<edges.size(); ++i)
	{
		real x = edges[i]->X;
		real y = edges[i]->Y;
		A[3*i] = x;
		A[3*i+1] = y;
		A[3*i+2] = 1.0;
		b[i] = -(x*x + y*y);
	}
	vector<real> z = LeastSquaresFitting(A, b, 3, bFlag);
	if(bFlag)
	{
		real sum = 0;
		for(int i=0; i<edges.size(); ++i)
		{
			real x = edges[i]->X;
			real y = edges[i]->Y;
			real err = x*x + y*y + z[0]*x + z[1]*y + z[2];
			sum += err*err;
		}
		return sum / ((z[0]*z[0]+z[1]*z[1])/4.0-z[2]); //normalize by the radius squared
	}
	else
	{
		return 0;
	}
}

int 
numActives(const vector<EdgelQueue>& elist)
{
	int count = 0;
	for(int i=0; i<elist.size(); ++i)
	{
		Edgel* pe = elist[i].Top();
		if(pe && pe->Active)
		{
			count++;
		}
	}
	return count;
}

vector<Edgel*>
collectActives(const vector<EdgelQueue>& elist)
{
	vector<Edgel*> actives;
	for(int i=0; i<elist.size(); ++i)
	{
		Edgel* pe = elist[i].Top();
		if(pe && pe->Active)
		{
			if(find(actives.begin(), actives.end(), pe) == actives.end())
			{
				actives.push_back(pe);
			}
		}
	}
	return actives;
}

real
computeFittingError(const vector<Edgel*>& edges,
					const CParticleF& click,
					bool& bFlag)
{
	if(!enclosingClick(edges, click))
	{
		bFlag = false;
		return 0;
	}
	//real error = computeCircularFitError(edges, bFlag) / edges.size();
	real error = computeElastica(edges, bFlag) / edges.size();
	if(bFlag)
	{
		return error;
	}
	else
	{
		return 0;
	}
}

vector<real> 
computeDistances(const vector<EdgelQueue>&elist, 
				 const CParticleF& click, 
				 const int* dims)
{
	real maxD = sqrt((real)dims[0]*dims[0] + dims[1]*dims[1]);
	vector<real> vdata(elist.size());
	for(int i=0; i<elist.size(); ++i)
	{
		Edgel* pe = elist[i].Top();
		if(pe == NULL)
		{
			vdata[i] = maxD;
		}
		else
		{
			vdata[i] = Distance(click.m_X, click.m_Y, pe->X, pe->Y);
		}
	}
	return vdata;
}

#include <szIndexedData.h>

vector<int> 
indexByProximity(const vector<real>& vdistances) 
{
	vector<indexedData> vdata(vdistances.size());
	for(int i=0; i<vdistances.size(); ++i)
	{
		vdata[i].index = i;
		vdata[i].data = vdistances[i]; 
	}
	sort(vdata.begin(), vdata.end());

	vector<int> vindex(vdata.size());
	for(int i=0; i<vindex.size(); ++i)
	{
		vindex[i] = vdata[i].index;
	}
	return vindex;
}

void
updateOrder(vector<int>& vindex,
			vector<real>& vdistance,
			int k)
{
	real distance = vdistance[k];
	if(k==0)
	{
		while(k<vdistance.size()-1 && vdistance[k]>vdistance[k+1])
		{
			swap(vdistance[k], vdistance[k+1]);
			swap(vindex[k], vindex[k+1]);
			k++;
		}
	}
	else if(k == vdistance.size() - 1)
	{
		while(k>0 && vdistance[k]<vdistance[k-1])
		{
			swap(vdistance[k], vdistance[k-1]);
			swap(vindex[k], vindex[k-1]);
			k--;
		}
	}
	else if(vdistance[k] < vdistance[k-1])
	{
		do
		{
			swap(vdistance[k], vdistance[k-1]);
			swap(vindex[k], vindex[k-1]);
			k--;
		}
		while(k>0 && vdistance[k]<vdistance[k-1]);
	}
	else if(vdistance[k] > vdistance[k+1])
	{
		do
		{
			swap(vdistance[k], vdistance[k+1]);
			swap(vindex[k], vindex[k+1]);
			k++;
		}
		while(k<vdistance.size()-1 && vdistance[k]>vdistance[k+1]);
	}
}

