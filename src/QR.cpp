#include <QR.h>
#include <cmath>
#include <szMexUtility.h>
#include <mex.h>

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*
Taken from Numerical Recipe.
*/
void qrdcmp(float **a, int n, float *c, float *d, int& sing)
{
	int i,j,k;
	float scale,sigma,sum,tau;

	sing=0;
	for (k=1;k<n;k++) {
		scale=0.0;
		for (i=k;i<=n;i++) scale=Max(scale,fabs(a[i][k]));
		if (scale == 0.0) {
			sing=1;
			c[k]=d[k]=0.0;
		} else {
			for (i=k;i<=n;i++) a[i][k] /= scale;
			for (sum=0.0,i=k;i<=n;i++) sum += SQR(a[i][k]);
			sigma=SIGN(sqrt(sum),a[k][k]);
			a[k][k] += sigma;
			c[k]=sigma*a[k][k];
			d[k] = -scale*sigma;
			for (j=k+1;j<=n;j++) {
				for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j];
				tau=sum/c[k];
				for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
			}
		}
	}
	d[n]=a[n][n];
	if (d[n] == 0.0) sing=1;
}

#undef SQR
#undef SIGN

std::vector<real> _Householder(const std::vector<real>& u, real c, int n)
{
	std::vector<real> M(n*n, 0);
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			M[i*n+j] = (i==j ? 1.0: 0) - u[i]*u[j]/c; 
		}
	}
	return M;
}

std::vector<real> _MultiplySquareMatrices(const std::vector<real>& a, const std::vector<real>& b, int n)
{
	std::vector<real> c(n*n, 0);
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			real sum = 0;
			for(int k=0; k<n; ++k)
			{
				sum += a[i*n+k] * b[k*n+j];
			}
			c[i*n+j] = sum;
		}
	}
	return c;
}

bool QR(const std::vector<real>& A, std::vector<real>& Q, std::vector<real>& R, int n)
{
	real** a = new real*[n+1];
	real* c = new real[n+1];
	real* d = new real[n+1];
	for(int i=0; i<=n; ++i)
	{
		a[i] = new real[n+1];
	}
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			a[i+1][j+1] = A[i*n+j];
		}
	}
	int sing;
	qrdcmp(a, n, c, d, sing);

	Q.resize(n*n, 0);
	R.resize(n*n, 0);

	for(int i=0; i<n; ++i)
	{
		for(int j=i; j<n; ++j)
		{
			R[i*n+j] = a[i+1][j+1];
		}
		R[i*n+i] = d[i+1];
		Q[i*n+i] = 1.0;
	}
	for(int j=0; j<n-1; ++j)
	{
		std::vector<real> u(n, 0);
		for(int i=j; i<n; ++i)
		{
			u[i] = a[i+1][j+1];
		}
		Q = _MultiplySquareMatrices(_Householder(u, c[j+1], n), Q, n);
	}
	//I have to transpose Q... I don't know why...
	for(int i=0; i<n; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			swap(Q[i*n+j], Q[j*n+i]);
		}
	}

	delete [] c;
	delete [] d;
	for(int i=0; i<=n; ++i)
	{
		delete [] a[i];
	}
	delete [] a;

	return sing ? false: true;
}

