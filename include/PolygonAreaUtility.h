#pragma once
#include <szParticleF.h>
#include <vector>
using namespace std;

vector<vector<CParticleF>>
triangulateInsidePolygon(vector<CParticleF> shape, bool& bSuccess);

float
polygonArea(vector<vector<CParticleF>>& triangles);

vector<vector<CParticleF>>
polygonOverlap(vector<vector<CParticleF>>& polygon1, vector<vector<CParticleF>>& polygon2);

bool
boundingBox(vector<CParticleF>& polygon, float& minx, float& miny, float& maxx, float& maxy);