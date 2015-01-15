#include <szContour.h>
#include <mexFileIO.h>

int
LoadContour(vector<Contour>& contours, const mxArray *prhs)
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
			int ndim2 = mxGetNumberOfDimensions(ar);
			const int* dims2 = (const int*) mxGetDimensions(ar);
			float* pf = (float*) mxGetData(ar);
			Contour c;
			c.points.resize(dims2[0]);
			for(int j=0; j<dims2[0]; ++j)
			{
				c.points[j].m_X = pf[j];
				c.points[j].m_Y = pf[j+dims2[0]];
			}
			contours.push_back(c);
		}
		return 0;
	}
	else
	{
		mexPrintf("LoadContour(): contours must be in cell array.\n");
		return MyMexUnsupportedClassError; //return it anyway
	}
}

mxArray*
StoreContours(const vector<Contour>& contours)
{
	const int dims[] = {contours.size()};
	mxArray* cell = mxCreateCellArray(1, (const mwSize*) dims);
	for(int i=0; i<contours.size(); ++i)
	{
		int n = contours[i].points.size();
		const int dimsC[] = {n, 3};
		mxArray* ar = mxCreateNumericArray(2, (const mwSize*) dimsC, mxSINGLE_CLASS, mxREAL);
		float* p = (float*) mxGetData(ar);
		for(int j=0; j<n; ++j)
		{
			p[j] = (float)contours[i].points[j].m_X;
			p[n+j] = (float)contours[i].points[j].m_Y;
			p[2*n+j] = (float)contours[i].points[j].m_Life;
		}
		mxSetCell(cell, i, ar);
	}
	return cell;
}
