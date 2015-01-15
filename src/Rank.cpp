//#include <Rank.h>

template<class T>
int
_partition(std::vector<T>& ar, int a, int b)
{
	T pivot = ar[b];
	int i=a;
	for(int j=a; j<b; ++j)
	{
		if(ar[j] < pivot)
		{
			swap(ar[i], ar[j]);
			i++;
		}
	}
	swap(ar[b], ar[i]);
	return i;
}

/*
The function find the location of the rank-k number in the 
vector, v.  The contents of v will be rearranged.
To find the median call the function with k=N/2 (or N/2 and (N-1)/2) 
if N is even) where N is the size of the vector.
*/
template<class T>
T
_Rank(std::vector<T>& ar, int i, int j, int k)
{
	int m = _partition(ar, i, j);
	if(m == k)
	{
		return ar[m];
	}
	else if(k < m)
	{
		return _Rank(ar, i, m-1, k);
	}
	else
	{
		return _Rank(ar, m+1, j, k);
	}
}

/*
A wrapper for the above Rank.
It performs the Rank filter on all data points.
*/
template<class T>
T
Rank(const std::vector<T>& ar, int order)
{
	std::vector<T> b = ar;
	return _Rank(b, 0, b.size()-1, order);
}
