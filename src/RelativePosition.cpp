#include <RelativePosition.h>
#include <cmath>
#include <szMexUtility.h>

RelativePosition
Side(float orientation, float direction)
{
	if(orientation != orientation || direction != direction)
	{
		return UnknownSide;
	}
	else
	{
		float df = direction - orientation;
		float sn = sin(df);
		float cs = cos(df);
		if(cs > sqrt(2.0)/2.0) return FrontSide;
		else if(cs < -sqrt(2.0)/2.0) return BackSide;
		else if(sn > 0) return LeftSide;
		else return RightSide;
	}
}

RelativePosition
Side(ContourEQW& a, ContourEQW& b)
{
	RelativePosition rp = UndeterminedSide;
	for(int i=0; i<a.size(); ++i)
	{
		float th = a.Orientation(i);
		float x = a.X[i];
		float y = a.Y[i];
		for(int j=0; j<b.size(); ++j)
		{
			float x2 = b.X[j];
			float y2 = b.Y[j];
			float th2 = atan2(y2-y, x2-x);
			RelativePosition rp2 = Side(th, th2);
			if(rp == UndeterminedSide)
			{
				rp = rp2;
			}
			else if(rp != rp2)
			{
				return UndeterminedSide;
			}
		}
	}
	return rp;
}

RelativePosition
FrontBack(ContourEQW& a, ContourEQW& b)
{
	RelativePosition rp = UndeterminedSide;
	float eps = 1.0e-6;
	for(int i=0; i<a.size(); ++i)
	{
		float th = a.Orientation(i);
		float x = a.X[i];
		float y = a.Y[i];
		for(int j=0; j<b.size(); ++j)
		{
			float x2 = b.X[j];
			float y2 = b.Y[j];
			if(Abs(y2-y)<eps && Abs(x2-x)<eps) continue;

			float th2 = atan2(y2-y, x2-x);
			RelativePosition rp2 = Side(th, th2);
			if(rp == UndeterminedSide)
			{
				rp = rp2;
				continue;
			}
			else if(rp2 == FrontSide || rp2 == BackSide)
			{
				if(rp != rp2)
				{
					return UndeterminedSide;
				}
			}
			else
			{
				return UndeterminedSide; //neither front nor back
			}
		}
	}
	return rp;
}
