#ifndef ___SZ_EDGEL_H___
#define ___SZ_EDGEL_H___
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;
#include <myDataType.h>

enum EdgelCompareEnum {EdgelStrength, EdgelFitness, EdgelAttribute};

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
		if(compareType == EdgelStrength)
		{
			return Strength < p.Strength;
		}
		else if(compareType == EdgelFitness)
		{
			return Fitness < p.Fitness;
		}
		else if(compareType == EdgelAttribute)
		{
			return Attribute < p.Attribute;
		}
	}
	bool operator <= (const EdgelEQW& p) const
	{
		if(compareType == EdgelStrength)
		{
			return Strength <= p.Strength;
		}
		else if(compareType == EdgelFitness)
		{
			return Fitness <= p.Fitness;
		}
		else if(compareType == EdgelAttribute)
		{
			return Attribute <= p.Attribute;
		}
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
		printf("%d %d %f %f %f %d\n", X, Y, Strength, Fitness, Attribute, State);
	}
	int X;
	int Y;
	int Z;
	int State;
	real Strength;
	real Angle;
	real Fitness;
	real Attribute;  //some general purpose storage
	bool Active;
	static EdgelCompareEnum compareType;
	friend struct EdgelFactory;

private:
	EdgelEQW(int x=0, int y=0, int z=0, real s=0, real f=0, real a=0, bool v=false)
	{
		X = x;
		Y = y;
		Z = z;
		Strength = s;
		Fitness = f;
		Attribute = a;
		Active = v;
		State = 0;
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

	EdgelEQW* Make(int x, int y, int z, real s=0, real f=0, real a=0)
	{
		EdgelEQW* p = new EdgelEQW(x, y, z, s, f, a);
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



#endif /* ___SZ_EDGEL_H___ */