#ifndef ___MY_EIGEN_H___
#define ___MY_EIGEN_H___

#include <vector>
using namespace std;

/*bool
Eigen(const vector<double>& C, 
	  vector<double>& eigenVal, 
	  vector<double>& eigenVec, 
	  int ndim);*/

/*
Compute eigenvalues and eigenvectors of a square matrix M.
*/

struct struct_eigen
{
	struct_eigen(int ndim=0)
	{
		if(ndim>0)
		{
			evals = vector<double>(ndim, 0);
			evecs = vector<vector<double>>(ndim, evals);
		}
	}
	vector<double> evals;
	vector<vector<double>> evecs;
};

struct_eigen
computeEigenValues(const vector<double>& M, int ndim);

#endif /* ___MY_EIGEN_H___ */