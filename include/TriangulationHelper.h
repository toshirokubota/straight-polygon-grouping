#ifndef ___TRIANGULATION_HELPER_H___
#define ___TRIANGULATION_HELPER_H___
#include <vector>
using namespace std;
#include <Triangulation.h>
#include <szParticleF.h>
#include <ContourEQW.h>
#include <FragmentInfo.h>

std::vector<CParticleF> 
vector2particle(const std::vector<float>& P, 
				const int* dims);


std::vector<Triangulation::_Internal::_indexed_triangle> 
indices2structs(const std::vector<CParticleF>& P, 
				const std::vector<int>& T, 
				const int* dimsT);

std::vector<std::vector<FragmentInfo>>
TriangulationPointsFixedInterval(const std::vector<ContourEQW>& contours, int interval);

std::vector<std::vector<FragmentInfo>>
TriangulationPointsSampleSize(const std::vector<ContourEQW>& contours, int nump);

float
flatness(Triangulation::_Internal::_triangle* tr, Triangulation::_Internal::_edge* e);

float EffectiveLength(Triangulation::_Internal::_edge* e);

vector<CParticleF>
traceShape(const vector<Triangulation::_Internal::_triangle*> faces);

bool
isExposed(Triangulation::_Internal::_vertex* v);

Triangulation::_Internal::_vertex*
oppositeSideVertex(Triangulation::_Internal::_edge* e, Triangulation::_Internal::_triangle* tr);

Triangulation::_Internal::_vertex* 
oppositeSideVertexOfBoundaryEdge(Triangulation::_Internal::_edge* e);

Triangulation::_Internal::_vertex*
commonVertex(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2);

Triangulation::_Internal::_edge*
findEdge(Triangulation::_Internal::_vertex* u, Triangulation::_Internal::_vertex* v);

Triangulation::_Internal::_edge*
commonEdge(Triangulation::_Internal::_triangle* t, Triangulation::_Internal::_triangle* s);

float 
angleBetween(Triangulation::_Internal::_edge* e1, Triangulation::_Internal::_edge* e2);

void
RemoveDanglingEdges(Triangulation::Triangulator& trmap);

#endif /* ___TRIANGULATION_HELPER_H___ */