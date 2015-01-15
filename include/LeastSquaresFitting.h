#ifndef ___LEAST_SQUARE_FIT_H___
#define ___LEAST_SQUARE_FIT_H___
#include <vector>
#include <myDataType.h>
#include <GaussJordan.h>

/*
Perform the least square fitting of (Ax-b)^2 where A is a n by m matrix, 
b is a n-vector and x is a m-vector.
To populate A, A[0..numparams-1] is the first row, A[numparams, ..., 2*numparams-1]
is the second row, and so on.
*/
std::vector<double>
LeastSquaresFitting(const std::vector<double>& A, const std::vector<double>& b, int m, bool& bSuccess);

/*
Compute L2-norm of Ax-b.
*/
double
LeastSquaresFittingError(const std::vector<double>& A, const std::vector<double>& x, const std::vector<double>& b);

#endif /* ___LEAST_SQUARE_FIT_H___ */