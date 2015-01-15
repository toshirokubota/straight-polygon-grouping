#include <CircleFitting.h>
#include <LevenbergMarquardt.h>
#include <derivify.h>

double evalFunc(const vector<double>& vmodel, int n,  const void* param)
{
	//vmodel is [radius, centerx, centery]
	//params are points on a circumference.
	double* points = (double*) param;
	double x = points[2*n];
	double y = points[2*n+1];
	return vmodel[0]*vmodel[0] - ((x-vmodel[1])*(x-vmodel[1]) + (y-vmodel[2])*(y-vmodel[2]));
}

surreal evalFuncD(const vector<surreal>& vmodel, int n,  const void* param)
{
	//vmodel is [radius, centerx, centery]
	//params are points on a circumference.
	double* points = (double*) param;
	double x = points[2*n];
	double y = points[2*n+1];
	return vmodel[0]*vmodel[0] - ((x-vmodel[1])*(x-vmodel[1]) + (y-vmodel[2])*(y-vmodel[2]));
}

bool fitCircle(const vector<CParticleF>& points,
				double& radius,
				double& centerx, 
				double & centery,
				int maxIter,
				double rate)
{
	double* parray = new double[points.size()*2];
	vector<double> vvalues;
	for(int i=0; i<points.size(); ++i)
	{
		vvalues.push_back(0.0);
		parray[i*2] = (double)points[i].m_X;
		parray[i*2+1] = (double)points[i].m_Y;
	}
	vector<double> vmodel(3, 0.0);
	vmodel[0] = radius;
	vmodel[1] = centerx;
	vmodel[2] = centery;
	bool bSuccess = LevenbergMarquardt(vvalues, vmodel, evalFunc, evalFuncD, maxIter, rate, parray);

	delete [] parray;
	radius = abs(vmodel[0]);
	centerx = vmodel[1];
	centery = vmodel[2];
	return bSuccess;
}

double fittingCircleError(const vector<CParticleF>& points,
				double radius,
				double centerx, 
				double centery)
{
	double* parray = new double[points.size()*2];
	vector<double> vvalues;
	for(int i=0; i<points.size(); ++i)
	{
		vvalues.push_back(0.0);
		parray[i*2] = (double)points[i].m_X;
		parray[i*2+1] = (double)points[i].m_Y;
	}
	vector<double> vmodel(3, 0.0);
	vmodel[0] = radius;
	vmodel[1] = centerx;
	vmodel[2] = centery;
	double error = 0;
	for(int i=0; i<points.size(); ++i)
	{
		double er = evalFunc(vmodel, i, parray);
		error += er*er;
	}

	delete [] parray;
	return error;
}