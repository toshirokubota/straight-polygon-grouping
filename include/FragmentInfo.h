#ifndef ___FRAGMENT_INFO_H___
#define ___FRAGMENT_INFO_H___

/*
A light weight structure to carry a contour fragment and its orientation.
*/
struct FragmentInfo
{
	FragmentInfo(int id=0, int pi=0): contourID(id), pointID(pi), saliency(0), angle(0) {}
	int contourID;
	//bool orientation;
	int pointID;
	float saliency;
	float angle;
	bool operator==(const FragmentInfo& f) const
	{
		return f.pointID == pointID && f.contourID == contourID;
	}
};

inline float EdgeWeight(float d) 
{
	return d*d;
}

#endif /* ___FRAGMENT_INFO_H___ */