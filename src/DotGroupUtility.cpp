#include <DotGroupUtility.h>
#include <szDistanceTransform.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>

void
updateDistanceMap(vector<float>& dmap, 
				  const vector<CParticleF>& shape, 
				  const int* dims)
{
	vector<unsigned char> outer(dims[0]*dims[1], (unsigned char)1);
	for(int y=0; y<dims[1]; ++y)
	{
		for(int x=0; x<dims[0]; ++x)
		{
			if(inside(CParticleF((float)x, (float)y), shape))
			{
				SetData2(outer, x, y, dims[0], dims[1], (unsigned char) 0);
			}
		}
	}
	for(int i=0; i<shape.size(); ++i)
	{
		SetData2(outer, shape[i].m_X, shape[i].m_Y, dims[0], dims[1], (unsigned char) 0);
	}
	DistanceTransformF(dmap, outer, DistanceTransformMode_Euclid, 2, dims);
}
