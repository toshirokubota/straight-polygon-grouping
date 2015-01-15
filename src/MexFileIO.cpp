
//#include "mexFileIO.h"

template<class Item>
int
LoadData(vector<Item>& A, const mxArray *prhs, mxClassID& class_id, int& ndim, const int** dims) 
{
	ndim      = mxGetNumberOfDimensions(prhs);
	*dims      = (const int*) mxGetDimensions(prhs);
	int numelm = mxGetNumberOfElements(prhs);
	class_id = mxGetClassID(prhs);
	A = vector<Item> (numelm);
	if (A.size()<numelm) 
	{
		mexPrintf("Failed to allocate memory in LoadData().\n");
		return MyMexAllocationError; //return it anyway
	}

	if (class_id == mxCELL_CLASS)
	{
		mexPrintf("use LoadCellData for cell array.\n");
		A.clear();
		return MyMexUnsupportedClassError; //return it anyway
	}
	else if (class_id == mxDOUBLE_CLASS)
	{
		double* V = (double*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxSINGLE_CLASS)
	{
		float* V = (float*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxINT32_CLASS)
	{
		int* V = (int*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxUINT32_CLASS)
	{
		unsigned int* V = (unsigned int*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxINT16_CLASS)
	{
		short* V = (short*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxUINT16_CLASS)
	{
		unsigned short* V = (unsigned short*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxINT8_CLASS)
	{
		char* V = (char*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxCHAR_CLASS)
	{
		mxChar* V = mxGetChars(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxUINT8_CLASS)
	{
		unsigned char* V = (unsigned char*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else if (class_id == mxLOGICAL_CLASS)
	{
		bool* V = (bool*) mxGetData(prhs);	
		for(int i=0; i<numelm; ++i)
			A[i]=(Item) V[i];
	}
	else {
		mexPrintf("unsupported input class in LoadData().\n");
		mexPrintf("supported are double, single, int32, uint32, int16, uint16, int8, uint8 or logical.\n");
		mexPrintf("use LoadCellData for cell array.\n");
		A.clear();
		return MyMexUnsupportedClassError; //return it anyway
	}
	return 0;
}

template<class Item>
int
ReadScalar(Item& A, const mxArray *prhs, mxClassID& class_id) 
{
	int ndim      = mxGetNumberOfDimensions(prhs);
	const int* dims = (const int*) mxGetDimensions(prhs);
	int numelm = mxGetNumberOfElements(prhs);
	class_id = mxGetClassID(prhs);

	if(numberOfElements(ndim,dims)<1)
		return -1;

	if (class_id == mxDOUBLE_CLASS)
	{
		double* V = (double*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxSINGLE_CLASS)
	{
		float* V = (float*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxINT32_CLASS)
	{
		int* V = (int*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxUINT32_CLASS)
	{
		unsigned int* V = (unsigned int*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxINT16_CLASS)
	{
		short* V = (short*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxUINT16_CLASS)
	{
		unsigned short* V = (unsigned short*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxINT8_CLASS)
	{
		char* V = (char*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxUINT8_CLASS)
	{
		unsigned char* V = (unsigned char*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else if (class_id == mxLOGICAL_CLASS)
	{
		bool* V = (bool*) mxGetData(prhs);	
		A =(Item) V[0];
	}
	else 
	{
		mexPrintf("unsupported input class in LoadData().\n");
		mexPrintf("use double, single, int32, uint32, int16, uint16, int8, uint8 or logical\n.");
		return MyMexUnsupportedClassError; //return it anyway
	}
	return 0;
}

template<class Item>
mxArray*
StoreData(vector<Item>& vdata, mxClassID class_id, int ndim, const int* dims) 
{
	mxArray* plhs;
	plhs   = mxCreateNumericArray(ndim, (const mwSize*) dims, class_id, mxREAL);
	int numelm=numberOfElements(ndim,dims);
	//mexPrintf("StoreData: Number of elements = %d\n", numelm);

	if (class_id == mxDOUBLE_CLASS)
	{
		double* p2 = (double*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(double) vdata[i];
	}
	else if (class_id == mxSINGLE_CLASS)
	{
		float* p2 = (float*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(float) vdata[i];
	}
	else if (class_id == mxINT32_CLASS)
	{
		int* p2 = (int*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i) {
			p2[i]=(int) vdata[i];
		}
	}
	else if (class_id == mxUINT32_CLASS)
	{
		unsigned int* p2 = (unsigned int*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(unsigned int) vdata[i];
	}
	else if (class_id == mxINT16_CLASS)
	{
		short* p2 = (short*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(short) vdata[i];
	}
	else if (class_id == mxUINT16_CLASS)
	{
		unsigned short* p2 = (unsigned short*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(unsigned short) vdata[i];
	}
	else if (class_id == mxINT8_CLASS)
	{
		char* p2 = (char*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(char) vdata[i];
	}
	else if (class_id == mxUINT8_CLASS)
	{
		unsigned char* p2 = (unsigned char*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(unsigned char) vdata[i];
	}
	else if (class_id == mxLOGICAL_CLASS)
	{
		bool* p2 = (bool*) mxGetData(plhs);
		for(int i=0; i<numelm; ++i)
			p2[i]=(bool) vdata[i];
	}
	else {
		mexPrintf("unsupported input class in StoreData().\n");
		mexPrintf("use double, single, int32, uint32, int16, uint16, int8, uint8 or logical\n.");
		//return MyMexUnsupportedClassError; //return it anyway
	}
	return plhs;
}

template<class Item>
mxArray*
StoreScalar(Item data, mxClassID class_id) 
{
	mxArray* plhs;
	const int dims[1]={1};
	plhs   = mxCreateNumericArray(1, (const mwSize *)dims, class_id, mxREAL);
	//mexPrintf("StoreData: Number of elements = %d\n", numelm);

	if (class_id == mxDOUBLE_CLASS)
	{
		double* p2 = (double*) mxGetData(plhs);
		p2[0]=(double) data;
	}
	else if (class_id == mxSINGLE_CLASS)
	{
		float* p2 = (float*) mxGetData(plhs);
		p2[0]=(float) data;
	}
	else if (class_id == mxINT32_CLASS)
	{
		int* p2 = (int*) mxGetData(plhs);
		p2[0]=(int) data;
	}
	else if (class_id == mxUINT32_CLASS)
	{
		unsigned int* p2 = (unsigned int*) mxGetData(plhs);
		p2[0]=(unsigned int) data;
	}
	else if (class_id == mxINT16_CLASS)
	{
		short* p2 = (short*) mxGetData(plhs);
		p2[0]=(short) data;
	}
	else if (class_id == mxUINT16_CLASS)
	{
		unsigned short* p2 = (unsigned short*) mxGetData(plhs);
		p2[0]=(unsigned short) data;
	}
	else if (class_id == mxINT8_CLASS)
	{
		unsigned char* p2 = (unsigned char*) mxGetData(plhs);
		p2[0]=(unsigned char) data;
	}
	else if (class_id == mxUINT8_CLASS)
	{
		char* p2 = (char*) mxGetData(plhs);
		p2[0]=(char) data;
	}
	else if (class_id == mxLOGICAL_CLASS)
	{
		bool* p2 = (bool*) mxGetData(plhs);
		p2[0]=(bool) data;
	}
	else 
	{
		mexPrintf("unsupported input class in StoreData().\n");
		mexPrintf("use double, single, int32, uint32, int16, uint16, int8, uint8 or logical\n.");
		//return MyMexUnsupportedClassError; //return it anyway
	}
	return plhs;
}

/*
This function stores a vector of vector of numbers into an array of cells.
Use NCOLS to set the number of columns for each cell.
*/
template<class Item>
mxArray*
StoreDataCell(vector<vector<Item>>& vdata, mxClassID class_id, int ndim, const int* dims, int ncols) 
{
	mxArray* plhs;
	plhs   = mxCreateCellArray((mwSize) ndim, (const mwSize*) dims);
	int numelm=numberOfElements(ndim,dims);
	for(int i=0; i<numelm; ++i)
	{
		int dims0[2] = { vdata[i].size() / ncols, ncols};
		mxArray* a = StoreData(vdata[i], class_id, 2, dims0);
		mxSetCell(plhs, i, a);
	}
	return plhs;
}

