#include <iostream>
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
#include <set>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <Triangulation.h>
#include <FragmentInfo.h>
#include <TriangulationHelper.h>
#include <szDistanceTransform.h>
#include <DistanceMapUtility.h>
#include <DotGroupUtility.h>
#include <intersectionConvexPolygons.h>
#include <szConvexHull2D.h>
#include <TriangulationSaliencyUtility.h>
#include <Eigen.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	printf("%s: This build was compiled at %s %s\n", "Testbed", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [J] = Testbed(M)");
		return;
	}

	vector<double> M;
	mxClassID classId;
	int ndim;
	const int* dims;
	LoadData(M, prhs[0], classId, ndim, &dims);

	if(ndim != 2 || dims[0] != dims[1])
	{
		mexErrMsgTxt("The input is not a square matrix.\n");
	}
	struct_eigen st = computeEigenValues(M, dims[0]);
	printf("Eigenvalues:\n");
	for(int i=0; i<dims[0]; ++i)
	{
		printf("%f ", st.evals[i]);
	}
	printf("\n");
	printf("Eigenvectors:\n");
	for(int i=0; i<dims[0]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			printf("%f ", st.evecs[i][j]);
		}
		printf("\n");
	}

	mexUnlock();
}

