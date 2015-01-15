#ifndef ___SZ_EDGEL_EQW_H___
#define ___SZ_EDGEL_EQW_H___
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;
#include <myDataType.h>

//enum EdgelCompareEnum {EdgelStrength, EdgelFitness, EdgelAttribute};

struct EdgelEQW
{
public:
	bool operator == (const EdgelEQW& p) const
	{
		if(X == p.X && Y == p.Y && Z == p.Z)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator < (const EdgelEQW& p) const
	{
		return Strength < p.Strength;
	}
	bool operator <= (const EdgelEQW& p) const
	{
		return Strength <= p.Strength;
	}
	bool operator > (const EdgelEQW& p) const
	{
		if(p < *this)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator >= (const EdgelEQW& p) const
	{
		if(p<=*this)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	void Print() const
	{
		printf("%f %f %f %f %f %f\n", X, Y, A, B, C, D);
	}
	real X;
	real Y;
	real A;
	real B;
	real C;
	real D;
	real Strength;
	friend struct EdgelFactory;

private:
	EdgelEQW(real x=0, real y=0, real a=0, real b=0, real c=0, real d=0)
	{
		X = x;
		Y = y;
		A = a;
		B = b;
		C = c;
		D = d;
	}
};

struct EdgelQueue
{
	void Enqueue(EdgelEQW* p)
	{
		if(find(edges.begin(), edges.end(), p) == edges.end())
		{
			edges.push_back(p);
		}
	}
	void Dequeue()
	{
		if(edges.empty())
		{
			return;
		}
		else
		{
			edges.erase(edges.begin());
			return;
		}
	}
	EdgelEQW* Top() const {return edges.empty() ? NULL: edges[0];}
	int length() const {return edges.size();}
	bool Swap()
	{
		if(edges.size() <= 1)
		{
			return false;
		}
		else
		{
			swap(edges[0], edges[1]);
			return true;
		}
	}

	vector<EdgelEQW*> edges;
};

struct EdgelFactory
{
	static EdgelFactory* GetInstance()
	{
		if(_instance == 0)
		{
			_instance = new EdgelFactory();
		}
		return _instance;
	}

	EdgelEQW* Make(real x, real y, real a=0, real b=0, real c=0, real d=0)
	{
		EdgelEQW* p = new EdgelEQW(x, y, a, b, c, d);
		edges.push_back(p);
		return p;
	}
	bool Remove(EdgelEQW* pe)
	{
		vector<EdgelEQW*>::iterator p = find(edges.begin(), edges.end(), pe);
		if(p == edges.end())
		{
			return false;
		}
		else
		{
			edges.erase(p);
			delete pe;
			return true;
		}
	}
protected:
	EdgelFactory()
	{
	}
	~EdgelFactory()
	{
		for(int i=0; i<edges.size(); ++i)
		{
			delete edges[i];
		}
	}
private:
	static EdgelFactory* _instance;
	vector<EdgelEQW*> edges;
};



#endif /* ___SZ_EDGEL_EQW_H___ */