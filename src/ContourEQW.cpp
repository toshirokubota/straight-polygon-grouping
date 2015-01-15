#include <ContourEQW.h>
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
int
LoadContourEQW(vector<ContourEQW>& contours, const mxArray *prhs)
{
	int ndim = mxGetNumberOfDimensions(prhs);
	const int* dims = (const int*) mxGetDimensions(prhs);
	mxClassID class_id = mxGetClassID(prhs);
	int numelm = mxGetNumberOfElements(prhs);
	if (class_id == mxCELL_CLASS)
	{
		for(int i=0; i<numelm; ++i)
		{
			mxArray* ar = mxGetCell(prhs, i);

			int ndimA;
			const int* dimsA;
			mxClassID classA;
			vector<float> A;
			LoadData(A,ar,classA,ndimA,&dimsA);
			int n = dimsA[0];
			ContourEQW c(n);
			for(int j=0; j<n; ++j)
			{
				c.X[j] = GetData2(A, j, 0, dimsA[0], dimsA[1], (float)0);
				c.Y[j] = GetData2(A, j, 1, dimsA[0], dimsA[1], (float)0);
				c.A[j] = GetData2(A, j, 2, dimsA[0], dimsA[1], (float)0);
				c.B[j] = GetData2(A, j, 3, dimsA[0], dimsA[1], (float)0);
				c.C[j] = GetData2(A, j, 4, dimsA[0], dimsA[1], (float)0);
				c.D[j] = GetData2(A, j, 5, dimsA[0], dimsA[1], (float)0);
				c.Strength[j] = GetData2(A, j, 6, dimsA[0], dimsA[1], (float)0);
			}
			contours.push_back(c);
		}
		return 0;
	}
	else
	{
		mexPrintf("LoadContourEQW(): contours must be in cell array.\n");
		return MyMexUnsupportedClassError; //return it anyway
	}
}

mxArray*
StoreContoursEQW(const vector<ContourEQW>& contours)
{
	const int dims[] = {contours.size()};
	mxArray* cell = mxCreateCellArray(1, (mwSize*) dims);
	for(int i=0; i<contours.size(); ++i)
	{
		int n = contours[i].X.size();
		const int dimsC[] = {n, 7};
		mxArray* ar = mxCreateNumericArray(2, (mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = contours[i].X[j];
			p[n+j] = contours[i].Y[j];
			p[2*n+j] = contours[i].A[j];
			p[3*n+j] = contours[i].B[j];
			p[4*n+j] = contours[i].C[j];
			p[5*n+j] = contours[i].D[j];
			if(j < contours[i].Strength.size())
			{
				p[6*n+j] = contours[i].Strength[j];
			}
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}

ContourEQW::ContourEQW(int n):X(n,0), Y(n,0), A(n,0), B(n,0), C(n,0), D(n,0), Strength(n,0)
{
}

ContourEQW::ContourEQW(Contour c)
{
	int n = c.points.size();
	X.resize(n);
	Y.resize(n);
	A.resize(n);
	B.resize(n);
	C.resize(n);
	D.resize(n);
	Strength.resize(n);
	for(int i=0; i<c.points.size(); ++i)
	{
		X[i] = c.points[i].m_X;
		Y[i] = c.points[i].m_Y;
		A[i] = 0;
		B[i] = 0;
		C[i] = 0;
		D[i] = 0;
		Strength[i] = c.points[i].m_Life;
	}
}

ContourEQW
ContourEQW::operator =(Contour c)
{
	int n = c.points.size();
	X.resize(n);
	Y.resize(n);
	A.resize(n);
	B.resize(n);
	C.resize(n);
	D.resize(n);
	Strength.resize(n);
	for(int i=0; i<n; ++i)
	{
		X[i] = c.points[i].m_X;
		Y[i] = c.points[i].m_Y;
		A[i] = 0;
		B[i] = 0;
		C[i] = 0;
		D[i] = 0;
		Strength[i] = c.points[i].m_Life;
	}
	return *this;
}

vector<float>
ContourEQW:: _FindOptimumPoint(int n0, int n2) const 
{
	float x0[3]={X[n0], A[n0], B[n0]};
	float x2[3]={X[n2], A[n2], B[n2]};
	float v[3],x1[3];
	int i,j;
	for(i=0; i<3; ++i) 
	{
		float sum=.0;
		for(j=0; j<3; ++j) 
		{
			sum+=H[i][j]*x2[j]+H[j][i]*x0[j];
		}
		v[i]=sum;
	}
	for(i=0; i<3; ++i) 
	{
		float sum=.0;
		for(j=0; j<3; ++j) 
		{
			sum+=G[i][j]*v[j];
		}
		x1[i]=-sum;
	}
	vector<float> res(6);
	res[0]=x1[0];
	res[2]=x1[1];
	res[3]=x1[2];

	float y0[3]={Y[n0], C[n0], D[n0]};
	float y2[3]={Y[n2], C[n2], D[n2]};
	float u[3],y1[3];
	for(i=0; i<3; ++i) 
	{
		float sum=.0;
		for(j=0; j<3; ++j) 
		{
			sum+=H[i][j]*y2[j]+H[j][i]*y0[j];
		}
		u[i]=sum;
	}
	for(i=0; i<3; ++i) 
	{
		float sum=.0;
		for(j=0; j<3; ++j) 
		{
			sum+=G[i][j]*u[j];
		}
		y1[i]=-sum;
	}
	res[1]=y1[0];
	res[4]=y1[1];
	res[5]=y1[2];

	return res;	
}

void
ContourEQW::Smoothen(int niter, float r, float t, bool bClosed, bool bFixedEnds)
{
	float a=4*r+t*t, b=8*r+t*t, c=2*r+t*t;
	G[0][0] = .5*a/b; G[0][1]=-.5/b; G[0][2] = 0;
	G[1][0] = -.5/b; G[1][1] = 1./(t*t*b); G[1][2] = 0;
	G[2][0] = 0; G[2][1] = 0; G[2][2] = .5/c;
	H[0][0] = -2.; H[0][1]=-t*t; H[0][2] = t;
	H[1][0] = -t*t; H[1][1] = 0; H[1][2] = -2.*r*t;
	H[2][0] = -t; H[2][1] = 2.*r*t; H[2][2] = -2.*r;

	int nump=X.size();
	if(nump <= 1)
	{
		return; //needs at least 2 points.
	}
	for(int i=0; i<niter; ++i)
	{
		for(int j=0; j<nump; ++j)
		{
			vector<float>res;
			if(bClosed)
			{
				int j0 = (j == 0) ? nump - 1: j - 1;
				int j2 = (j == nump-1) ? 0: j + 1;
				res = _FindOptimumPoint(j0, j2);
			}
			else 
			{
				if(j==0)
				{
					if(bFixedEnds)
					{
					res.resize(6);
					res[0]=X[j]; 
					res[1]=Y[j];
					res[2]=0; 
					res[3]=X[j+1]-X[j];
					res[4]=0;
					res[5]=Y[j+1]-Y[j];
					}
					else
					{
						int j0 = 0;
						int j2 = (j == nump-1) ? 0: j + 1;
						res = _FindOptimumPoint(j0, j2);
					}
				}
				else if(j==nump-1)
				{
					if(bFixedEnds)
					{
					res.resize(6);
					res[0]=X[j]; 
					res[1]=Y[j];
					res[2]=0; 
					res[3]=X[j]-X[j-1];
					res[4]=0;
					res[5]=Y[j]-Y[j-1];
					}
					else
					{
						int j0 = (j == 0) ? 0: j - 1;
						int j2 = nump-1;
						res = _FindOptimumPoint(j0, j2);
					}
				}
				else
				{
					int j0 = j - 1;
					int j2 = j + 1;
					res = _FindOptimumPoint(j0, j2);
				}
			}
			X[j] = res[0];
			Y[j] = res[1];
			A[j] = res[2];
			B[j] = res[3];
			C[j] = res[4];
			D[j] = res[5];
		}
	}
}


float
ContourEQW::EvaluateX(int i, float t) const
{
	return X[i] + A[i]*t*t + B[i]*t;
}

float
ContourEQW::EvaluateY(int i, float t) const
{
	return Y[i] + C[i]*t*t + D[i]*t;
}

float
ContourEQW::Curvature(int i) const
{
	float dd = pow(B[i]*B[i]+D[i]*D[i], (float)1.5);
	return (B[i]*C[i] - A[i]*D[i])/Max(0.00001,dd);
}

float
ContourEQW::Orientation(int i) const
{
	//return atan2(D[i], B[i]);
	//return atan2(B[i], D[i]);
	if(i==0) return atan2(Y[i+1]-Y[i], X[i+1]-X[i]);
	else if(i==size()-1) return atan2(Y[i]-Y[i-1], X[i]-X[i-1]);
	else return atan2(Y[i+1]-Y[i-1], X[i+1]-X[i-1]);
}

vector<float>
ContourEQW::Curvatures() const
{
	vector<float> curvatures(X.size());
	for(int i=0; i<X.size(); ++i)
	{
		float dd = pow(B[i]*B[i]+D[i]*D[i], (float)1.5);
		curvatures[i] = (B[i]*C[i] - A[i]*D[i])/Max(0.00001,dd);
	}
	return curvatures;
}

ContourEQW
ContourEQW::extract(int i, int j) const
{
	int first = Max(0, Min(i,j));
	int last = Min(X.size()-1, Max(i,j));
	ContourEQW s(last-first+1);
	for(int i=first; i<=last; ++i)
	{
		s.A[i-first]=A[i];
		s.B[i-first]=B[i];
		s.C[i-first]=C[i];
		s.D[i-first]=D[i];
		s.X[i-first]=X[i];
		s.Y[i-first]=Y[i];
		s.Strength[i-first]=Strength[i];
	}
	return s;
}

ContourEQW
ContourEQW::reverse() const
{
	ContourEQW s(size());
	for(int i=size()-1, j=0; i>=0; --i, ++j)
	{
		s.A[j]=A[i];
		s.B[j]=-B[i]; //reverse the sign. I am not sure if this is correct...
		s.C[j]=C[i];
		s.D[j]=-D[i]; //reverse the sign. I am not sure if this is correct...
		s.X[j]=X[i];
		s.Y[j]=Y[i];
		s.Strength[j]=Strength[i];
	}
	return s;
}

ContourEQW
extractContourFragment(ContourEQW& c,
					   int first, int last)
{
	ContourEQW res;
	for(int i=first; i<=last; ++i)
	{
		res.A.push_back(c.A[i]);
		res.B.push_back(c.B[i]);
		res.C.push_back(c.C[i]);
		res.D.push_back(c.D[i]);
		res.X.push_back(c.X[i]);
		res.Y.push_back(c.Y[i]);
		res.Strength.push_back(c.Strength[i]);
	}
	return res;
}

