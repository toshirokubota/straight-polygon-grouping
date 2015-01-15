//#include <mex.h>
#include <vector>
//#include <algorithm>
//#include <queue>
#include <limits>
#include <map>
using namespace std;
//#include <szMyNeighborOp.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <szParticleF.h>
#include <Triangulation.h>
#include <TriangulationHelper.h>
//#include <DistanceMapUtility.h>
//#include <szConvexHull2D.h>

/*float
similarityMeasure(Triangulation::_Internal::_triangle* s, Triangulation::_Internal::_triangle* t)
{
	if(s==0 || t==0) return 0;

	Triangulation::_Internal::_edge* common = commonEdge(s, t);
	if(common == NULL) return 0;
	Triangulation::_Internal::_vertex* sv = oppositeSideVertex(common, s);
	Triangulation::_Internal::_edge* e1 = findEdge(sv, common->vertices[0]);
	Triangulation::_Internal::_edge* e2 = findEdge(sv, common->vertices[1]);
	float df = (common->Length() - Max(e1->Length(), e2->Length()));
	return 1./(1. + df*df);
}*/

/*float
similarityMeasure(Triangulation::_Internal::_triangle* s, Triangulation::_Internal::_triangle* t)
{
	if(s==0 || t==0) return 0;

	Triangulation::_Internal::_edge* common = commonEdge(s, t);
	if(common == NULL) return 0;
	Triangulation::_Internal::_vertex* sv = oppositeSideVertex(common, s);
	Triangulation::_Internal::_vertex* tv = oppositeSideVertex(common, t);
	float sa = areaTriangle(s->vertices[0]->p, s->vertices[1]->p, s->vertices[2]->p);
	float ta = areaTriangle(t->vertices[0]->p, t->vertices[1]->p, t->vertices[2]->p);
	float dif = Abs(sa - ta) / common->Length();

	vector<CParticleF> pnts;
	pnts.push_back(sv->p);
	pnts.push_back(tv->p);
	pnts.push_back(common->vertices[0]->p);
	pnts.push_back(common->vertices[1]->p);
	vector<CParticleF> hull = ConvexHull2D(pnts);
	float ha = polygonArea(hull);
	float dif2 = Max(0, ha - sa - ta);

	float eps = 100;
	return eps/(eps + dif + dif2);
}*.

/*float
similarityMeasure(Triangulation::_Internal::_triangle* s, Triangulation::_Internal::_triangle* t)
{
	if(s==0 || t==0) return 0;

	Triangulation::_Internal::_edge* common = commonEdge(s, t);
	if(common == NULL) return 0;
	Triangulation::_Internal::_edge* se[2];
	for(int i=0, j=0; i<3; ++i)
	{
		if(s->edges[i]!=common)
		{
			se[j++] = s->edges[i];
		}
	}
	Triangulation::_Internal::_edge* te[2];
	for(int i=0, j=0; i<3; ++i)
	{
		if(t->edges[i]!=common)
		{
			te[j++] = t->edges[i];
		}
	}
	float df1 = (Max(se[0]->Length(), se[1]->Length()) - Max(te[0]->Length(), te[1]->Length())) / common->Length();
	float df2 = (Min(se[0]->Length(), se[1]->Length()) - Min(te[0]->Length(), te[1]->Length())) / common->Length();
	return 1.0/(1.0 + df1*df1 + df2*df2);
}*/

float
similarityMeasure1(Triangulation::_Internal::_triangle* s, Triangulation::_Internal::_triangle* t)
{
	if(s==0 || t==0) return 0;

	Triangulation::_Internal::_edge* common = commonEdge(s, t);
	if(common == NULL) return 0;
	Triangulation::_Internal::_vertex* sv = oppositeSideVertex(common, s);
	Triangulation::_Internal::_vertex* tv = oppositeSideVertex(common, t);
	CParticleF pnts[4];
	pnts[0] = sv->p;
	pnts[1] = tv->p;
	pnts[2] = common->vertices[0]->p;
	pnts[3] = common->vertices[1]->p;
	CParticleF c;
	for(int i=0; i<4; ++i)
	{
		c.m_X += pnts[i].m_X;
		c.m_Y += pnts[i].m_Y;
	}
	c.m_X /= 4.0;
	c.m_Y /= 4.0;
	float d = Distance2LineSegment(pnts[2], pnts[3], c) / common->Length();
	return 1.0 / (1.0 + d);
}

float
similarityMeasure(Triangulation::_Internal::_triangle* s, Triangulation::_Internal::_triangle* t)
{
	if(s==0 || t==0) return 0;

	Triangulation::_Internal::_edge* common = commonEdge(s, t);
	float slen = Max(s->edges[0]->Length(), Max(s->edges[1]->Length(), s->edges[2]->Length()));
	float tlen = Max(t->edges[0]->Length(), Max(t->edges[1]->Length(), t->edges[2]->Length()));
	return common->Length() / Max(slen, tlen);
}

/*
This version does not assume two triangles are adjacent.
*/
float
similarityMeasure0(Triangulation::_Internal::_triangle* s, Triangulation::_Internal::_triangle* t)
{
	if(s==0 || t==0) return 0;
	vector<float> vsides1;
	for(int i=0; i<3; ++i)
	{
		vsides1.push_back(s->edges[i]->Length());
	}
	sort(vsides1.begin(), vsides1.end());
	vector<float> vsides2;
	for(int i=0; i<3; ++i)
	{
		vsides2.push_back(t->edges[i]->Length());
	}
	sort(vsides2.begin(), vsides2.end());
	float maxdif = Max(Abs(vsides1[0]-vsides2[0]), Max(Abs(vsides1[1]-vsides2[1]), Abs(vsides1[2]-vsides2[2])));
	float minlen = Min(vsides1[0], vsides2[0]);
	return minlen / (minlen + maxdif);
}

vector<float>
computeSimilarityMeasures(Triangulation::Triangulator& trmap)
{
	vector<float> measures(trmap.edges.size(), 0);
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_triangle* s = trmap.edges[i]->faces[0];
		Triangulation::_Internal::_triangle* t = trmap.edges[i]->faces[1];
		if(s && t)
		{
			measures[i] = similarityMeasure(s, t);
		}
	}
	return measures;
}

vector<int>
clusterBySimilarity(Triangulation::Triangulator& trmap, float thres)
{
	vector<Node<Triangulation::_Internal::_triangle*>*> nodes;
	map<Triangulation::_Internal::_triangle*,int> fmap;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		nodes.push_back(makeset(trmap.faces[i]));
		fmap[trmap.faces[i]] = i;
	}
	vector<float> measures = computeSimilarityMeasures(trmap);
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		Triangulation::_Internal::_triangle* s = trmap.edges[i]->faces[0];
		Triangulation::_Internal::_triangle* t = trmap.edges[i]->faces[1];
		if(s && t)
		{
			if(measures[i] > thres)
			{
				int k1 = fmap[s];
				int k2 = fmap[t];
				merge(nodes[k1], nodes[k2]);
			}
		}
	}
	vector<Node<Triangulation::_Internal::_triangle*>*> cls = clusters(nodes);
	vector<int> labels(trmap.faces.size(), 0);
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		labels[i] = distance(cls.begin(), find(cls.begin(), cls.end(), findset(nodes[i]))) + 1;
	}
	for(int i=0; i<nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}

vector<int>
clusterBySimilarityV2(Triangulation::Triangulator& trmap, 
						float thres, //threshold for adjacent triangle similarity
						float thres2 //threshold for non-adjacent triangle similarity
						)
{
	vector<Node<Triangulation::_Internal::_triangle*>*> nodes;
	map<Triangulation::_Internal::_triangle*,int> fmap;
	vector<vector<Node<Triangulation::_Internal::_triangle*>*>> groups;
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		Node<Triangulation::_Internal::_triangle*>* n = makeset(trmap.faces[i]);
		nodes.push_back(n);
		fmap[trmap.faces[i]] = i;
		vector<Node<Triangulation::_Internal::_triangle*>*> v(1, n);
		groups.push_back(v);
	}
	vector<float> measures = computeSimilarityMeasures(trmap);
	vector<pair<float,Triangulation::_Internal::_edge*>> pairs;
	for(int i=0; i<trmap.edges.size(); ++i)
	{
		pair<float,Triangulation::_Internal::_edge*> p(-measures[i], trmap.edges[i]);
		pairs.push_back(p);
	}
	sort(pairs.begin(), pairs.end());

	for(int i=0; i<trmap.edges.size(); ++i)
	{
		float measure = -pairs[i].first;
		Triangulation::_Internal::_edge* ed = pairs[i].second;
		Triangulation::_Internal::_triangle* s = ed->faces[0];
		Triangulation::_Internal::_triangle* t = ed->faces[1];
		if(s && t)
		{
			if(measure > thres)
			{
				int k1 = fmap[s];
				int k2 = fmap[t];
				Node<Triangulation::_Internal::_triangle*>* n1 = findset(nodes[k1]);
				Node<Triangulation::_Internal::_triangle*>* n2 = findset(nodes[k2]);
				if(n1 == n2) continue;

				int j1 = fmap[n1->key];
				int j2 = fmap[n2->key];
				//float minmes = std::numeric_limits<float>::infinity();
				int count = 0;
				for(int g1=0; g1<groups[j1].size(); ++g1)
				{
					Triangulation::_Internal::_triangle* tr1 = groups[j1][g1]->key;
					for(int g2=0; g2<groups[j2].size(); ++g2)
					{
						Triangulation::_Internal::_triangle* tr2 = groups[j2][g2]->key;
						float m = similarityMeasure0(tr1, tr2);
						if(m > thres2) count++;
					}
				}
				if(count > (groups[j1].size() * groups[j2].size())/2)
				{
					merge(nodes[k1], nodes[k2]);
					Node<Triangulation::_Internal::_triangle*>* n3 = findset(nodes[k1]);
					if(n3==n1)
					{
						groups[j1].insert(groups[j1].end(), groups[j2].begin(), groups[j2].end());
					}
					else
					{
						groups[j2].insert(groups[j2].end(), groups[j1].begin(), groups[j1].end());
					}
				}
			}
		}
	}

	vector<Node<Triangulation::_Internal::_triangle*>*> cls = clusters(nodes);
	vector<int> labels(trmap.faces.size(), 0);
	for(int i=0; i<trmap.faces.size(); ++i)
	{
		labels[i] = distance(cls.begin(), find(cls.begin(), cls.end(), findset(nodes[i]))) + 1;
	}
	for(int i=0; i<nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}
