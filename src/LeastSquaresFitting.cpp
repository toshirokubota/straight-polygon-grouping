#include <LeastSquaresFitting.h>

/*
Perform the least square fitting of (Ax-b)^2 where A is a n by m matrix, 
b is a n-vector and x is a m-vector.
*/
std::vector<double>
LeastSquaresFitting(const std::vector<double>& A, const std::vector<double>& b, int m, bool& bSuccess)
{
	std::vector<double> M(m*m);
	std::vector<double> y(m);
	int n = b.size();
	for(int i=0; i<m; ++i)
	{
		for(int j=0; j<m; ++j)
		{
			double sum = 0;
			for(int k=0; k<n; ++k)
			{
				sum += A[i+m*k] * A[j+m*k];
			}
			M[i*m+j] = sum;
		}
	}
	for(int i=0; i<m; ++i)
	{
		double sum = 0;
		for(int k=0; k<n; ++k)
		{
			sum += A[i+m*k] * b[k];
		}
		y[i] = sum;
	}

	bSuccess = GaussJordan(M, y, m);
	return y;
}

/*
Compute L2-norm of Ax-b.
*/
double
LeastSquaresFittingError(const std::vector<double>& A, const std::vector<double>& x, const std::vector<double>& b)
{
	double error = 0;
	int m = x.size();
	for(int i=0; i<b.size(); ++i)
	{
		double val = 0;
		for(int j=0; j<x.size(); ++j)
		{
			val += A[m*i+j] * x[j];
		}
		error += (val - b[i]) * (val - b[i]);
	}
	return error;
}