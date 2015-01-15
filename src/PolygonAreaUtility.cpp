#include <PolygonAreaUtility.h>
#include <szMiscOperations.h>
#include <IntersectionConvexPolygons.h>

namespace _PolygonAreaUtility
{
	vector<CParticleF>
	reverseShape(vector<CParticleF>& shape)
	{
		vector<CParticleF>  rev(shape.size());
		for (int i = 0, j = 0; i < shape.size(); ++i, j = j == 0 ? shape.size() - 1 : j - 1)
		{
			rev[i] = shape[j];
		}
		return rev;
	}

	bool
	isEar(int idx, vector<CParticleF>& shape)
	{
		CParticleF a = shape[idx == 0 ? shape.size() - 1 : idx - 1];
		CParticleF b = shape[idx];
		CParticleF c = shape[idx == shape.size() - 1 ? 0 : idx + 1];
		CParticleF m((a.m_X + c.m_X) / 2, (a.m_Y + c.m_Y) / 2);
		if (inside(m, shape) == false)
		{
			return false;
		}
		for (int j = 0; j < shape.size(); ++j)
		{
			if (j == idx || j == (idx - 2 + shape.size()) % shape.size() || j == (idx - 1 + shape.size()) % shape.size() || j == (idx + 1) % shape.size())
				continue;
			CParticleF p = shape[j];
			CParticleF q = shape[(j + 1) % shape.size()];
			pair<float, float> para = _IntersectConvexPolygon::intersect(a, c, p, q);
			if (para.first > 0 && para.first < 1.0 && para.second > 0 && para.second < 1.0)
			{
				return false;
			}
		}
		return true;
	}
}

vector<vector<CParticleF>>
triangulateInsidePolygon(vector<CParticleF> shape, bool& bSuccess)
{
	if (ClockWise(shape) < 0)
	{
		shape = _PolygonAreaUtility::reverseShape(shape);
	}
	bSuccess = true;
	vector<vector<CParticleF>> triangles;
	while (shape.size() >= 3)
	{
		bool bFound = false;
		for (int i = 0; i < shape.size(); ++i)
		{
			if (shape.size() == 3 || _PolygonAreaUtility::isEar(i, shape))
			{
				int idx = i;
				CParticleF a = shape[idx == 0 ? shape.size() - 1 : idx - 1];
				CParticleF b = shape[idx];
				CParticleF c = shape[idx == shape.size() - 1 ? 0 : idx + 1];
				vector<CParticleF> tr(3);
				tr[0] = a; tr[1] = b; tr[2] = c;
				triangles.push_back(tr);
				shape.erase(shape.begin() + idx);
				bFound = true;
				break;
			}
		}
		if (bFound == false)
		{
			printf("Incomplete exit:\n");
			for (int j = 0; j < shape.size(); ++j)
			{
				printf("%f %f\n", shape[j].m_X, shape[j].m_Y);
			}
			bSuccess = false;
			break;
		}
	}
	return triangles;
}

float
polygonArea(vector<vector<CParticleF>>& triangles)
{
	float area = 0, area2 = 0, area3 = 0;
	for (int i = 0; i < triangles.size(); ++i)
	{
		area += areaTriangle(triangles[i][0], triangles[i][1], triangles[i][2]);
	}
	return area;
}

vector<vector<CParticleF>>
polygonOverlap(vector<vector<CParticleF>>& polygon1, vector<vector<CParticleF>>& polygon2)
{
	vector<vector<CParticleF>> overlap;
	for (int i = 0; i < polygon1.size(); ++i)
	{
		for (int j = 0; j < polygon2.size(); ++j)
		{
			vector<CParticleF> cv = _IntersectConvexPolygon::intersectConvexHulls(polygon1[i], polygon2[j]);
			if (cv.size() > 2)
			{
				float a = polygonArea(cv);
				if (a > 0)
				{
					overlap.push_back(cv);
				}
			}
		}
	}
	return overlap;
}

bool
boundingBox(vector<CParticleF>& polygon, float& minx, float& miny, float& maxx, float& maxy)
{
	if (polygon.empty()) return false;
	minx = maxx = polygon[0].m_X;
	miny = maxy = polygon[0].m_Y;

	for (int i = 1; i < polygon.size(); ++i)
	{
		float x = polygon[i].m_X;
		float y = polygon[i].m_Y;
		if (minx > x) minx = x;
		if (miny > y) miny = y;
		if (maxx < x) maxx = x;
		if (maxy < y) maxy = y;
	}
	return true;
}

