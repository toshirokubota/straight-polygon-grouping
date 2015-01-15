#ifndef ___TRIANGULATION_SALIENCY_UTILITY_H___
#define ___TRIANGULATION_SALIENCY_UTILITY_H___
#include <Triangulation.h>
#include <szMexUtility.h>


float lengthTerm(float length, float lengthThres, float alpha, float angleTerm=0.0f);

float angleTerm(float angle);

float angleTerm(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2);

float 
pairFitness(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2, float beta, float gamma);

vector<Triangulation::_Internal::_vertex*>
orientVertices(vector<Triangulation::_Internal::_vertex*>& vertices, CParticleF& o);

/*
Given a treble of vertices (u, v, w) with u being in the middle, compute the saliency of u.
*/
float
pointFitness(Triangulation::_Internal::_vertex* u, 
			 Triangulation::_Internal::_vertex* v,
			 Triangulation::_Internal::_vertex* w,
			 float beta, float gamma);

vector<float>
Moment(vector<CParticleF>& pnts);

vector<float>
Moment2(vector<CParticleF>& pnts);

vector<float>
MomentCircular(vector<CParticleF>& pnts);

float
pointRadialFitness(Triangulation::_Internal::_vertex* u, 
					float beta, float gamma);

struct StructuralTensor
{
	StructuralTensor()
	{
		cx = 0; cy = 0;
		major = 0;
		minor = 0;
		major_axis[0] = 0; major_axis[1] = 0;
		minor_axis[0] = 0; minor_axis[1] = 0;
	}
	StructuralTensor(vector<CParticleF>& pnts)
	{
		vector<float> eg = Moment2(pnts);
		cx = eg[6]; cy = eg[7];
		if(Abs(eg[0]) > Abs(eg[1]))
		{
			major = Abs(eg[0]); minor=Abs(eg[1]);
			major_axis[0] = eg[2]; major_axis[1] = eg[3];
			minor_axis[0] = eg[4]; minor_axis[1] = eg[5];
		}
		else
		{
			major = Abs(eg[1]); minor=Abs(eg[0]);
			major_axis[0] = eg[4]; major_axis[1] = eg[5];
			minor_axis[0] = eg[2]; minor_axis[1] = eg[3];
		}
	}
	float cx, cy;
	float major, minor;
	float major_axis[2];
	float minor_axis[2];
};

struct StructuralTensor3D
{
	StructuralTensor3D()
	{
		_init();
	}
	StructuralTensor3D(vector<CParticleF>& pnts);
	void _init()
	{
		cx = 0; cy = 0; cz = 0;
		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
			{
				axes[i][j] = 0;
			}
			evals[i] = 0;
		}
	}
	float cx, cy, cz;
	float evals[3];
	float axes[3][3];
};

#endif /* ___TRIANGULATION_SALIENCY_UTILITY_H___ */
