#include <ContourEQWRegistration.h>
#include <szParticleF.h>
#include <LeastSquaresFitting.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMiscOperations.h>
#include <QR.h>

real
ArcLength(const ContourEQW& c, int beg, int end)
{
	real length = 0;
	for(int i=Max(0, beg); i<=end-1 && i<c.size()-1; ++i)
	{
		length += Distance(c.X[i], c.Y[i], c.X[i+1], c.Y[i+1]);
	}
	return length;
}

real
SimilarityByTranslation(const ContourEQW& a, const ContourEQW& b)
{
	CParticleF ca = Centroid(a);
	CParticleF cb = Centroid(b);
	real d = Distance(ca, cb);

	real value = d * d; // + d * (a.size() - b.size()) * (a.size() - b.size());

	return value;
}

real
SimilarityByLength(const ContourEQW& a, const ContourEQW& b)
{
	real value = ArcLength(a, 0, a.size()-1) - ArcLength(b, 0, b.size()-1);
	return Abs(value);
}

real
SimilarityByOrientation(const ContourEQW& a, const ContourEQW& b)
{
	real t1 = atan2(a.Y[0] - a.Y[a.size()-1], a.X[0] - a.X[a.size()-1]);
	real t2 = atan2(b.Y[0] - b.Y[b.size()-1], b.X[0] - b.X[b.size()-1]);
	return 1.0 - cos(t1-t2)*cos(t1-t2);
}

ContourEQW
MoveToOrigin(const ContourEQW& a)
{
	CParticleF ac = Centroid(a);
	ContourEQW b = a;
	for(int i=0; i<b.size(); ++i)
	{
		b.X[i] -= ac.m_X;
		b.Y[i] -= ac.m_Y;
	}
	return b;
}

real
HausdorffDistance(const ContourEQW& a, const ContourEQW& b)
{
	real dist = 0;
	for(int i=0; i<a.size(); ++i)
	{
		real d0 = std::numeric_limits<real>::infinity();
		for(int j=0; j<b.size(); ++j)
		{
			real d = Distance(a.X[i], a.Y[i], b.X[j], b.Y[j]);
			if(d < d0) d0 = d;
		}
		if(d0 > dist) dist = d0;
	}
	for(int i=0; i<b.size(); ++i)
	{
		real d0 = std::numeric_limits<real>::infinity();
		for(int j=0; j<a.size(); ++j)
		{
			real d = Distance(a.X[j], a.Y[j], b.X[i], b.Y[i]);
			if(d < d0) d0 = d;
		}
		if(d0 > dist) dist = d0;
	}
	return dist;
}

real
TotalDifference(const ContourEQW& a, const ContourEQW& b)
{
	real dist = 0;
	{
		real inc = (real)b.size()/(real)a.size();
		real s=0;
		for(int i=0; i<a.size(); ++i, s+=inc)
		{
			real d = Distance(a.X[i], a.Y[i], b.EvaluateX(floor(s), s-floor(s)), b.EvaluateY(floor(s), s-floor(s)));
			dist += d;
		}
	}
	{
		real inc = (real)a.size()/(real)b.size();
		real s=0;
		for(int i=0; i<b.size(); ++i, s+=inc)
		{
			real d = Distance(b.X[i], b.Y[i], a.EvaluateX(floor(s), s-floor(s)), a.EvaluateY(floor(s), s-floor(s)));
			dist += d;
		}
	}
	return dist;
}

/*
Use linear registration to align contour B to contour A.
*/
real
ComputeSimilarityOneDirection(const ContourEQW& a, const ContourEQW& b)
{
	vector<real> x;
	ContourEQW c = ProcrustesRegistration(a, b, x);
	real sum = 0; //HausdorffDistance(b, c);
	/*for(int i=0; i<b.size(); ++i)
	{
		sum += Distance(b.X[i], b.Y[i], c.X[i], c.Y[i]);
	}
	sum += TotalDifference(a, c);
	sum /= ArcLength(c, 0, c.size()-1);*/
	real al = ArcLength(a, 0, a.size()-1);
	real bl = ArcLength(b, 0, b.size()-1);
	real cl = ArcLength(c, 0, c.size()-1);
	//sum += (x[0]*x[0] + x[1]*x[1]);
	//sum += 10*acos(cos(2*x[4]))*bl;
	sum += acos(cos(2*x[4]));
	//sum += bl * Max(x[2]/x[3], x[3]/x[2]);
	//sum += (al-bl)*(al-bl);
	//real hd = HausdorffDistance(b, c);
	//sum += hd*hd;

	return sum; // / bl;
}

real
KullbeckLeibler(vector<real>& p, vector<real>& q)
{
	real sum = 0;
	for(int i=0; i<p.size(); ++i)
	{
		if(q[i]<=0 || p[i]<=0) continue;
		sum += p[i] * log(p[i]/q[i]);
	}
	return sum;
}

vector<int>
TransitionHistogram(const ContourEQW& a)
{
	vector<int> h(9, 0);
	for(int i=0; i<a.size()-1; ++i)
	{
		real dx = a.X[i+1]-a.X[i];
		real dy = a.Y[i+1]-a.Y[i];
		int idx = dx>0.5 ? 1: (dx<-0.5 ? -1: 0);
		int idy = dy>0.5 ? 1: (dy<-0.5 ? -1: 0);
		int k = (idy+1) * 3 + idx+1;
		h[k]++;
	}
	return h;
}

real
StructuralMeasure(const ContourEQW& a, const ContourEQW& b)
{
	real eps = 1.0e-3;
	vector<int> binsA = TransitionHistogram(a);
	vector<int> binsB = TransitionHistogram(b);
	vector<real> pdfA(binsA.size(), 0);
	{
		int count = 0;
		for(int i=0; i<binsA.size(); ++i)
		{
			count += binsA[i];
		}
		for(int i=0; i<pdfA.size(); ++i)
		{
			pdfA[i] = ((real)binsA[i]+eps) /((real)count + eps*binsA.size());
		}
	}
	vector<real> pdfB(binsB.size(), 0);
	{
		int count = 0;
		for(int i=0; i<binsB.size(); ++i)
		{
			count += binsB[i];
		}
		for(int i=0; i<pdfB.size(); ++i)
		{
			pdfB[i] = ((real)binsB[i]+eps) /((real)count + eps*binsB.size());
		}
	}
	return KullbeckLeibler(pdfA, pdfB);
}

real
ComputeSimilarity(const ContourEQW& a, const ContourEQW& b)
{
	return Min(StructuralMeasure(a, b), StructuralMeasure(a, b.reverse()));
	//return SimilarityByTranslation(a, b); //by offset
	//return SimilarityByLength(a, b);
	//return SimilarityByOrientation(a, b);
	//return ComputeSimilarityOneDirection(a, b);

	real error1 = ComputeSimilarityOneDirection(a, b);
	real error2 = ComputeSimilarityOneDirection(a, b.reverse());

	return Min(error1, error2);
}

