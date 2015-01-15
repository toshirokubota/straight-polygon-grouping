#include <ContourEQWRegistration.h>
#include <szmexutilitytemplate.h>
#include <LeastSquaresFitting.h>

ContourEQW
Interpolate(const ContourEQW& a)
{
	ContourEQW b(a.size()*2 - 1);
	for(int i=0; i<a.size()-1; ++i)
	{
		b.X[2*i] = a.X[i];
		b.Y[2*i] = a.Y[i];
		b.X[2*i+1] = (a.X[i]+a.X[i+1])/2.0;
		b.Y[2*i+1] = (a.Y[i]+a.Y[i+1])/2.0;
	}
	b.X[b.size()-1] = a.X[a.size()-1];
	b.Y[b.size()-1] = a.Y[a.size()-1];
	return b;
}


/*
Using a rigid transformation, register contour B to contour A.
The resulting contour has the same length with B. 
Only X and Y are set and coefficients are set to 0.
*/
ContourEQW
AffineRegistration(const ContourEQW& a, const ContourEQW& b, vector<real>& params)
{
	const int dims[] = {6, 2*b.size()};
	vector<real> v(dims[1]);
	vector<real> A(dims[0]*dims[1], 0);
	float inc = (float)a.size()/(float)b.size();
	float t = 0;
	for(int i=0; i<b.size(); ++i, t+=inc)
	{
		int k = (int)t;
		v[2*i] = a.EvaluateX(k, t-k);
		v[2*i+1] = a.EvaluateY(k, t-k);
		SetData2(A, 0, 2*i, dims[0], dims[1], b.X[i]);
		SetData2(A, 1, 2*i, dims[0], dims[1], b.Y[i]);
		SetData2(A, 2, 2*i, dims[0], dims[1], (float)1);
		SetData2(A, 3, 2*i+1, dims[0], dims[1], b.X[i]);
		SetData2(A, 4, 2*i+1, dims[0], dims[1], b.Y[i]);
		SetData2(A, 5, 2*i+1, dims[0], dims[1], (float)1);
	}
	bool bSuccess;
	vector<real> x = LeastSquaresFitting(A, v, 6, bSuccess);
	ContourEQW c = b;
	if(bSuccess)
	{
		for(int i=0; i<b.size(); ++i)
		{
			c.X[i] = x[0] * b.X[i] + x[1] * b.Y[i] + x[2];
			c.Y[i] = x[3] * b.X[i] + x[4] * b.Y[i] + x[5];
			c.A[i] = c.B[i]= c.D[i] = c.Strength[i] = 0;
		}
		params = x;
	}
	return c;
}

ContourEQW
AffineRegistration(const ContourEQW& a, const ContourEQW& b)
{
	vector<real> params;
	return AffineRegistration(a, b, params);
}

/*
Using a rigid transformation, register contour B to contour A.
The resulting contour has the same length with B. 
No translation is involved. No shear is involved.
Only X and Y are set and coefficients are set to 0.
*/
ContourEQW
EuclideanRegistration(const ContourEQW& a, const ContourEQW& b, vector<real>& params)
{
	const int dims[] = {4, 2*b.size()};
	vector<real> v(dims[1]);
	vector<real> A(dims[0]*dims[1], 0);
	float inc = (float)a.size()/(float)b.size();
	float t = 0;
	for(int i=0; i<b.size(); ++i, t+=inc)
	{
		int k = (int)t;
		v[2*i] = a.EvaluateX(k, t-k);
		v[2*i+1] = a.EvaluateY(k, t-k);
		SetData2(A, 0, 2*i, dims[0], dims[1], b.X[i]);
		SetData2(A, 1, 2*i, dims[0], dims[1], b.Y[i]);
		SetData2(A, 2, 2*i+1, dims[0], dims[1], b.X[i]);
		SetData2(A, 3, 2*i+1, dims[0], dims[1], b.Y[i]);
	}
	bool bSuccess;
	vector<real> x = LeastSquaresFitting(A, v, 4, bSuccess);
	ContourEQW c = b;
	if(bSuccess)
	{
		for(int i=0; i<b.size(); ++i)
		{
			c.X[i] = x[0] * b.X[i] + x[1] * b.Y[i];
			c.Y[i] = x[2] * b.X[i] + x[3] * b.Y[i];
			c.A[i] = c.B[i]= c.D[i] = c.Strength[i] = 0;
		}
		params = x;
	}
	return c;
}

ContourEQW
EuclideanRegistration(const ContourEQW& a, const ContourEQW& b)
{
	vector<real> params;
	return EuclideanRegistration(a, b, params);
}

/*
Find the centroid of the contour points.
*/
CParticleF
Centroid(const ContourEQW& c)
{
	CParticleF cent(0, 0);
	for(int i=0; i<c.size(); ++i)
	{
		cent.m_X += c.X[i];
		cent.m_Y += c.Y[i];
	}
	cent.m_X /= (real)c.size();
	cent.m_Y /= (real)c.size();
	return cent;
}

void
PCA(const ContourEQW& c, real& tx, real& ty, real& sxx, real& sxy, real& syy)
{
	CParticleF cent = Centroid(c);
	tx = cent.m_X;
	ty = cent.m_Y;
	sxx=0; sxy=0; syy=0;
	for(int i=0; i<c.size(); ++i)
	{
		real dx = c.X[i] - tx;
		real dy = c.Y[i] - ty;
		sxx += dx * dx;
		sxy += dx * dy;
		syy += dy * dy;
	}
	sxx = (sxx/(real)c.size());
	sxy = (sxy/(real)c.size());
	syy = (syy/(real)c.size());
	if(Abs(sxy) < 1.0e-5)
	{
		sxy = sxy < 0 ? -1.0e-5: 1.0e-5;
	}
}

/*
Using a Procrustes transformation, register contour B to contour A.
Transform consists of translation, uniform scaling, and rotation.
The resulting contour has the same length with B. 
Only X and Y are set and coefficients are set to 0.
*/
ContourEQW
ProcrustesRegistration(const ContourEQW& a, const ContourEQW& b, vector<real>& z)
{
	vector<real> pa(5);
	{
		real tx, ty, sxx, sxy, syy;
		PCA(a, tx, ty, sxx, sxy, syy);
		real det = (sxx-syy)*(sxx-syy) + 4*sxy*sxy;
		real eval1 = (sxx + syy + sqrt(det))/2.0;
		real eval2 = (sxx + syy - sqrt(det))/2.0;
		real theta = atan2((real)1.0, -(-sxx+syy-sqrt(det))/(2*sxy));
		pa[0] = tx; 
		pa[1] = ty;
		pa[2] = eval1 < 1.0e-5 ? 1: eval1;
		pa[3] = eval2 < 1.0e-5 ? 1: eval2;
		pa[4] = theta;
	}
	vector<real> pb(5);
	{
		real tx, ty, sxx, sxy, syy;
		PCA(b, tx, ty, sxx, sxy, syy);
		real det = (sxx-syy)*(sxx-syy) + 4*sxy*sxy;
		real eval1 = (sxx + syy + sqrt(det))/2.0;
		real eval2 = (sxx + syy - sqrt(det))/2.0;
		real theta = atan2((real)1.0, -(-sxx+syy-sqrt(det))/(2*sxy));
		pb[0] = tx; 
		pb[1] = ty;
		pb[2] = eval1 < 1.0e-5 ? 1: eval1; 
		pb[3] = eval2 < 1.0e-5 ? 1: eval2;
		pb[4] = theta;
	}
	z.resize(5);
	z[0] = pa[0] - pb[0];
	z[1] = pa[1] - pb[1];
	z[2] = sqrt(pa[2]/pb[2]);
	z[3] = sqrt(pa[3]/pb[3]);
	z[4] = pa[4] - pb[4];
	//if((int)b.Strength[0] == 59)
	/*{
		printf("%d: %3.3f, %3.3f, %3.3f, %3.3f, %3.3f:: %3.3f, %3.3f, %3.3f, %3.3f, %3.3f\n", 
			(int)b.Strength[0],
			pa[0], pa[1], pa[2], pa[3], pa[4], pb[0], pb[1], pb[2], pb[3], pb[4]);
	}*/

	ContourEQW c = b;
	for(int i=0; i<b.size(); ++i)
	{
		real x = b.X[i] - pb[0]; //shift it toward the origin
		real y = b.Y[i] - pb[1];
		real x2 = x*cos(pb[4])+y*sin(pb[4]); //rotate by -theta
		real y2 = -x*sin(pb[4])+y*cos(pb[4]);
		real x3 = x2 * z[2];	//scale it
		real y3 = y2 * z[3];
		real x4 = x3*cos(pa[4])-y3*sin(pa[4]); //rotate by theta of A
		real y4 = x3*sin(pa[4])+y3*cos(pa[4]);
		c.X[i] = x4 + pa[0]; //shift it toward A
		c.Y[i] = y4 + pa[1];

		//if you want to consider only translation, uncomment these
		//c.X[i] = x + pa[0]; //shift it toward A
		//c.Y[i] = y + pa[1];

		c.A[i] = c.B[i]= c.D[i] = c.Strength[i] = 0;
	}
	return c;
}

ContourEQW
ProcrustesRegistration(const ContourEQW& a, const ContourEQW& b)
{
	vector<real> params;
	return ProcrustesRegistration(a, b, params);
}

