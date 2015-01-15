#ifndef ___INDEXED_DATA_H____
#define ___INDEXED_DATA_H____

class indexedData  
{
public:
	indexedData(double d=0, int i=0)
	{
		data = d;
		index = i;
	}
	double data;
	int index;
};

inline bool 
operator<(const indexedData& s1, const indexedData& s2) {
	return s1.data < s2.data;
}

#endif /* ___INDEXED_DATA_H____ */