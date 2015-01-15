#include <szParticle.h>

/* 
*/
template<class Item>
void
LocalMaximum(vector<unsigned char>& M, 
			 const vector<Item>& V, 
			 const vector<unsigned char>& L, //ROI mask
			 const vector<int>& nbh,
			 bool strict,
			 int ndim,
			 const int* dims)
{

	int nvoxels = numberOfElements(ndim,dims);
	if(M.size()<nvoxels || V.size()<nvoxels || L.size()<nvoxels) 
	{
		//mexPrintf("An input image does not contain enough data.\n");
		return;
	}

	int i;
	//construct a subscript representation of the neighborhood
	//for efficient boundary check
	vector<int> vcoords;
	for(i=0; i<nbh.size(); ++i) 
	{
		vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
		vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
	}

	for(i=0; i<nvoxels; ++i) 
	{
		unsigned char lb=L[i];
		if(!lb)
		{
			SetData(M,i,(unsigned char) 0);
			continue; //not inside ROI
		}
		bool bLM=true;
		Item origV=V[i];

		vector<int> vsub = Ind2Sub(i,ndim,dims);
		for(int j=0; j<nbh.size(); ++j) 
		{
			if(NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims)) 
			{ 
				int k=i+nbh[j];
				if(L[k]) 
				{
					Item v2=V[k];
					if(strict) 
					{// strictly maximum
						if(v2>=origV) 
						{
							bLM=false;
							break;
						}
					}
					else 
					{
						if(v2>origV) 
						{
							bLM=false;
							break;
						}
					}
				}
			}
		}

		if(bLM) 
		{
			SetData(M,i,(unsigned char) 1);
		}
		else 
		{
			SetData(M,i,(unsigned char) 0);
		}
	}
}


/* 
*/
template<class Item>
void
LocalMinimum(vector<unsigned char>& M, 
			 const vector<Item>& V, 
			 const vector<unsigned char>& L,
			 const vector<int>& nbh,
			 bool strict,
			 int ndim,
			 const int* dims) 
{
	int nvoxels = numberOfElements(ndim,dims);
	if(M.size()<nvoxels || V.size()<nvoxels || L.size()<nvoxels) 
	{
		//mexPrintf("An input image does not contain enough data.\n");
		return;
	}

	int i;
	//construct a subscript representation of the neighborhood
	//for efficient boundary check
	vector<int> vcoords;
	for(i=0; i<nbh.size(); ++i) 
	{
		vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
		vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
	}

	for(i=0; i<nvoxels; ++i) 
	{
		unsigned char lb=L[i];
		if(!lb)
			continue; //not inside ROI
		bool bLM=true;
		bool b0;
		Item origV=V[i];

		vector<int> vsub = Ind2Sub(i,ndim,dims);
		for(int j=0; j<nbh.size(); ++j)
		{
			if(NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims))
			{ 
				int k=i+nbh[j];
				if(L[k])
				{
					bool b;
					Item v2=V[k];
					if(strict) {// strictly maximum
						if(v2<=origV)
						{
							bLM=false;
							break;
						}
					}
					else 
					{
						if(v2<origV) 
						{
							bLM=false;
							break;
						}
					}
				}
			}
		}

		if(bLM) 
		{
			SetData(M,i,(unsigned char) 1);
		}
		else 
		{
			SetData(M,i,(unsigned char) 0);
		}
	}
}

