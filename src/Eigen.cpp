#include <Eigen.h>

#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

namespace _Eigen
{
	/*
	Use Jacobi Transformations of a symetric matrix to obtain the 
	normalized eigenvectors and eigen values 
	-- Taken from NUMERICAL RECIPES IN C (pp 364-366) --
	*/

#define EPS 0.000000001
#define MAX_ITER 50

#define ROTATE(a, i, j, k, l, n) g=a[i*n+j]; h=a[k*n+l]; a[i*n+j] = g-s*(h+g*tau); \
	a[k*n+l] = h+s*(g-h*tau);

	/* 
	Computes all eigenvalues and eigenvectors of a real symmetric matrix.
	a: a vector made up by concatenating row vectors of the matrix.
	the content will be destroyed.
	n: the dimension of the symmetric matrix
	d: returns eigenvalues
	v: returns eigenvectors - arranged in a vector form
	nrot: the number of iterations executed in jacobi.
	*/

	void jacobi(double* a, int n, double* d, double* v, int* nrot) {

		int j,iq,ip,ip_times_n,i ;
		double tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*vector();

		b = new double[n];
		z = new double[n];


		for(ip_times_n=0, ip=0; ip<n; ++ip, ip_times_n+=n)  
		{

			/* Initialize the identity matrix */
			for(iq=0; iq<n; ++iq)v[ip_times_n + iq] = 0.0 ;
			v[ip_times_n + ip] = 1.0 ;

			/* Initailize b and d to be diagonal of a */
			b[ip] = d[ip] = a[ip_times_n + ip];
			z[ip] = 0.0 ;
		}

		*nrot = 0 ;
		for(i=0;i<MAX_ITER;++i)
		{
			/* Sum off-diagonal elements */
			sm=0.0 ;

			for(ip_times_n=0,ip=0;ip<n-1;ip++,ip_times_n+=n)
				for(iq=ip+1;iq<n;iq++)
					sm += fabs(a[ip_times_n + iq]);

			/*  If we have converged,  free the working vectors and return. */
			if(sm == 0.0)
			{
				delete [] b;
				delete [] z;
				return;
			}

			/* tresh is different on first three iterations...*/
			tresh=(i<3) ? 0.2*sm/(n*n) : 0.0 ;

			for(ip_times_n=0,ip=0;ip<n-1;ip++,ip_times_n+=n)
			{
				for(iq=ip+1;iq<n;++iq)
				{
					g=100.0*fabs(a[ip_times_n + iq]);

					/* After four sweeps, skip the rotation if the off-diagonal element is small */
					/* This test is taken directly from the text and looks a little suspect to me... */

					if(i > 3 && g < EPS)
						a[ip_times_n + iq] = 0.0 ;

					else if(fabs(a[ip_times_n+iq]) > tresh) 
					{
						h=d[iq]-d[ip];
						if(g < EPS)
							t = (fabs(a[ip_times_n+iq]) > EPS) ? (a[ip_times_n+iq])/h : 0.0 ; 
						else
						{ 
							theta=(fabs(h) < EPS) ? 0.0 : 0.5*h/(a[ip_times_n+iq]);
							t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
							if(theta < 0.0)
								t = -t ; 
						} 
						c=1.0/sqrt(1.0+t*t);
						s=t*c;
						tau=s/(1.0+c);

						h=t*a[ip_times_n+iq];
						z[ip] -= h;
						z[iq] += h;
						d[ip] -= h;
						d[iq] += h;
						a[ip_times_n+iq]=0.0;

						for(j=0;j<ip;j++)
						{
							ROTATE(a,j,ip,j,iq,n);
						}
						for(j=ip+1;j<iq;j++)
						{
							ROTATE(a,ip,j,j,iq,n);
						}
						for(j=iq+1;j<n;j++)
						{
							ROTATE(a,ip,j,iq,j,n);
						}
						for(j=0;j<n;j++)
						{
							ROTATE(v,j,ip,j,iq,n);
						}
						++(*nrot);
					}
				}
			}
			for(ip=0;ip<n;++ip)
			{
				b[ip] += z[ip];
				d[ip]=b[ip];
				z[ip]=0.0;
			}
		}

		/* Failed to converge!! What to do ??? */
		/* Well, let's at least free up memory and return without a murmur */


		delete [] b;
		delete [] z;
		return;	
	}

	/*
	A wrapper for jacobi
	*/
	bool
		Eigen(const vector<double>& C, 
		vector<double>& eigenVal, 
		vector<double>& eigenVec, 
		int ndim) 
	{

		int nrot;
		double* pC = new double[C.size()];
		double* pEval = new double[eigenVal.size()];
		double* pEvec = new double[eigenVec.size()];
		unsigned int i, j;
		for(i=0; i<C.size(); ++i)
		{
			pC[i] = C[i];
		}

		jacobi(pC, ndim, pEval, pEvec, &nrot);

		for(i=0; i<eigenVal.size(); ++i)
		{
			eigenVal[i] = (double)pEval[i];
		}
		for(i=0; i<ndim; ++i)
		{
			for(j=0; j<ndim; ++j)
			{
				eigenVec[i*ndim+j] = (double)pEvec[i+ndim*j];
			}
		}

		delete [] pC;
		delete [] pEval;
		delete [] pEvec;

		if(nrot == MAX_ITER)
			return false;
		else
			return true;
	}

}

struct_eigen
computeEigenValues(const vector<double>& M, int ndim)
{

	vector<double> evals(ndim);
	vector<double> evecs(ndim*ndim);

	_Eigen::Eigen(M, evals, evecs, ndim);
	vector<pair<double,int>> pairs;
	for(int i=0; i<evals.size(); ++i)
	{
		pair<double,int> p(-evals[i], i); //MINUS SIGN to sort in descending order
		pairs.push_back(p);
	}
	sort(pairs.begin(), pairs.end());
	struct_eigen st(ndim);
	for(int i=0; i<pairs.size(); ++i)
	{
		int k = pairs[i].second;
		st.evals[i] = evals[k];
		for(int j=0; j<ndim; ++j)
		{
			st.evecs[i][j] = evecs[k*ndim+j];
		}
	}
	return st;
}
