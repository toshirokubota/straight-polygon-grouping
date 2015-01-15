#include <stdio.h>
#include <cmath>

#include <vector>
#include <algorithm>
using namespace std;
#include <szMexUtility.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <TraceIterator.h>

struct DPEntry
{
	DPEntry(int t=0, int s=0, int val=0, DPEntry* p=NULL)
	{
		target = t;
		source = s;
		value = val;
		prev = p;
	}
	float value;
	DPEntry* prev;
	int target;
	int source;
};

vector<vector<DPEntry*>>
RunDP(vector<int>& t, vector<int>& s, float skip, float repeat, float (*Score)(int, int))
{
	vector<DPEntry*> row(t.size(), 0);
	vector<vector<DPEntry*>> table(s.size(), row);
	for(int i=0; i<table.size(); ++i)
	{
		for(int j=0; j<row.size(); ++j)
		{
			table[i][j] = new DPEntry(t[j], s[i]);
		}
	}

	table[0][0]->value = Score(s[0], t[0]);
	float sum = 0;
	for(int i=0; i<table.size(); ++i)
	{
		for(int j=0; j<row.size(); ++j)
		{
			if(i==0 && j==0) continue;

			float cost = std::numeric_limits<float>::infinity();
			float score = Score(s[i], t[j]);
			DPEntry* pi = NULL;
			//using a match of previous source symbol
			if(j>0 && i>0)
			{
				float val = table[i-1][j-1]->value + score;
				if(val < cost) 
				{
					cost = val;
					pi = table[i-1][j-1];
				}
			}
			if(i>0)
			{
				float val = table[i-1][j]->value + repeat + score;
				if(val < cost)
				{
					cost = val;
					pi = table[i-1][j];
				}
			}
			//using current source symbol repeating to match multiple target symbols
			if(j>0)
			{
				float val = table[i][j-1]->value + skip + score;
				if(val < cost)
				{
					cost = val;
					pi = table[i][j-1];
				}
			}
			table[i][j]->value = cost;
			table[i][j]->prev = pi;
		}
	}

	return table;
}

/*
*
* Random pattern measures
*
*/

/*
Given three points on a fragment, categorize the turn into one of five:
Front or stright, LeftFront or slight turn to left, RightFront or slight turn to right, 
LeftSide or sharp turn to left, and RightSide or sharp turn to right.
*/
int
getWalkDirection(int x0, int y0, int x1, int y1, int x2, int y2)
{
	TraceIterator it(x1, y1, x2, y2);
	int x, y;
	int k=0;
	while(it.Next(x, y))
	{
		if(Sign(x) == Sign(x0-x1) && Sign(y)==Sign(y0-y1))
		{
			break;
		}
		k++;
	}
	return k;
}

/*
Record a contour fragment on an integer grid. Skip those that remain in the same grid.
*/
vector<CParticleF>
gridTrace(const ContourEQW& c)
{
	vector<CParticleF> trace;
	if(c.size()==0) return trace;

	CParticleF prev((int)c.X[0], (int)c.Y[0]);
	trace.push_back(prev);
	for(int i=1; i<c.size(); ++i)
	{
		int x = (int)c.X[i];
		int y = (int)c.Y[i];
		if(x != (int)prev.m_X || y!=(int)prev.m_Y)
		{
			CParticleF curr(x, y);
			trace.push_back(curr);
			prev = curr;
		}
	}
	return trace;
}

vector<int>
Contour2WalkDirections(const ContourEQW& c)
{
	vector<CParticleF> trace = gridTrace(c);
	vector<int> index;
	if(trace.size()>2)
	{
		for(int j=2; j<trace.size(); ++j)
		{
			index.push_back(getWalkDirection(trace[j].m_X, trace[j].m_Y, trace[j-1].m_X, trace[j-1].m_Y, trace[j-2].m_X, trace[j-2].m_Y));
		}
	}
	return index;
}

float 
ScoreRandom(int t, int s)
{
	int map[]={1, 0, 2, 7, 3, 6, 4, 5}; //maps a trace-iterator index to a consecutive one in clock-wise starting from the top-left corner.
	return Min(Abs(map[t]-map[s]), 8-Abs(map[t]-map[s]));
}

float
RandomDesimilarityDP(const ContourEQW& tc, const ContourEQW& sc, float skip, float repeat)
{
	vector<int> t = Contour2WalkDirections(tc);
	vector<int> s = Contour2WalkDirections(sc);
	vector<vector<DPEntry*>> table = RunDP(t, s, skip, repeat, ScoreRandom);
	DPEntry* e = table[table.size()-1][table[table.size()-1].size()-1];
	float value = e->value;
	for(int i=0; i<table.size(); ++i)
	{
		for(int j=0; j<table[i].size(); ++j)
		{
			delete table[i][j];
		}
	}

	return value;
}

/*
*
* Structural patterns measure
*
*/

vector<int>
Contour2TransitionIndex(const ContourEQW& c)
{
	vector<int> v;
	float x1 = c.X[0];
	float y1 = c.Y[0];
	for(int i=1; i<c.size(); ++i)
	{
		float x2 = c.X[i];
		float y2 = c.Y[i];

		int iy=0, ix=0;
		if(y2-y1 < -.5) iy=-1;
		else if(y2-y1>0.5) iy = 1;
		if(x2-x1 < -.5) ix=-1;
		else if(x2-x1>0.5) ix = 1;

		if(iy==-1 && ix==-1) 
		{
			v.push_back(0);
		}
		else if(iy==-1 && ix==0)
		{
			v.push_back(1);
		}
		else if(iy==-1 && ix==1)
		{
			v.push_back(2);
		}
		else if(iy==0 && ix==-1) 
		{
			v.push_back(7);
		}
		else if(iy==0 && ix==0)
		{
			//do nothing
		}
		else if(iy==0 && ix==1)
		{
			v.push_back(3);
		}
		else if(iy==1 && ix==-1) 
		{
			v.push_back(6);
		}
		else if(iy==1 && ix==0)
		{
			v.push_back(5);
		}
		else if(iy==1 && ix==1)
		{
			v.push_back(4);
		}
		x1 = x2;
		y1 = y2;
	}
	return v;
}

float 
ScoreStructure(int t, int s)
{
	return Min(Abs(t-s), 8-Abs(t-s));
}

float
StructuralDesimilarityDP(const ContourEQW& tc, const ContourEQW& sc, float skip, float repeat)
{
	vector<int> t = Contour2TransitionIndex(tc);
	vector<int> s = Contour2TransitionIndex(sc);
	vector<vector<DPEntry*>> table = RunDP(t, s, skip, repeat, ScoreStructure);
	DPEntry* e = table[table.size()-1][table[table.size()-1].size()-1];
	float value = e->value;
	for(int i=0; i<table.size(); ++i)
	{
		for(int j=0; j<table[i].size(); ++j)
		{
			delete table[i][j];
		}
	}

	return value;
}
