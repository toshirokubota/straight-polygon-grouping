#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>

#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
#include <set>
#include <hash_map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szConvexHull2D.h>
#include <ContourEQW.h>
#include <szParticleF.h>
#include <DisjointSet.h>
#include <IntersectionConvexPolygons.h>
#include <GraphFactory.h>
#include <MovingParticle.h>
#include <ParticleSimulator.h>

GraphFactory<CParticleF>* GraphFactory<CParticleF>::_instance = NULL;
ParticleFactory* ParticleFactory::_instance = NULL;

int MovingParticle::_id = 0;

vector<pair<int,int>>
indices2pairs(vector<int> T, const int* dims)
{
	int n = dims[0];
	vector<pair<int,int>> pairs(n);
	for(int i=0; i<n; ++i)
	{
		int k1 = GetData2(T, i, 0, dims[0], dims[1], 0);
		int k2 = GetData2(T, i, 1, dims[0], dims[1], 0);
		pairs[i] = pair<int,int>(k1-1, k2-1);
	}
	return pairs;
}

vector<Snapshot>
Connect()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<float> evt;
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		(*it)->updateEvent();
		evt.push_back((*it)->getEvent().t);
	}
	sort(evt.begin(), evt.end());
	float thres = evt[0]*10;
	printf("thres = %f\n", thres);
	float time = 0;
	vector<Snapshot> parts;
	while (time < thres)
	{
		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;

		EventStruct ev = p->getEvent();
		if (ev.type == SplitEvent && ev.q != NULL && ev.r != NULL) {
			MovingParticle* q = (MovingParticle*)ev.q;
			MovingParticle* r = q->getNext();
			MovingParticle* pnb[2] = { p->getNext(), p->getPrev() };
			MovingParticle* qnb[2] = { q->getNext(), q->getPrev() };
			MovingParticle* rnb[2] = { r->getNext(), r->getPrev() };
			MovingParticle* p1 = factory->makeParticle(p->project(ev.t), Split, 0.0f);
			p1->setVelocity(0, 0);
			MovingParticle* p2 = factory->makeParticle(p->project(ev.t), Split, 0.0f);
			p2->setVelocity(0, 0);
			MovingParticle::setNeighbors(p1, p, r);
			MovingParticle::setNeighbors(p2, q, pnb[0]);
			vector<MovingParticle*> vp1 = MovingParticle::vectorize(p1);
			vector<MovingParticle*> vp2 = MovingParticle::vectorize(p2);

			vector<vector<MovingParticle*>> areas1 = MovingParticle::closedRegions(vp1);
			vector<vector<MovingParticle*>> areas2 = MovingParticle::closedRegions(vp2);
			for (int i = 0; i < areas1.size(); ++i) {
				parts.push_back(Snapshot(0.0f, areas1[i]));
				vector<CParticleF> pnts;
				for (int j = 0; j < areas1[i].size(); ++j) 
				{
					pnts.push_back(areas1[i][j]->getP0());
				}
				//if (ClockWise(pnts)==false)
				{
					for (int j = 0; j < areas1[i].size(); ++j) 
					{
						factory->inactivate(areas1[i][j]);
					}
				}
			}
			for (int i = 0; i < areas2.size(); ++i) {
				parts.push_back(Snapshot(0.0f, areas2[i]));
				vector<CParticleF> pnts;
				for (int j = 0; j < areas2[i].size(); ++j) 
				{
					pnts.push_back(areas2[i][j]->getP0());
				}
				//if (ClockWise(pnts)==false)
				{
					for (int j = 0; j < areas2[i].size(); ++j) 
					{
						factory->inactivate(areas2[i][j]);
					}
				}
			}

			//restore the state
			/*factory->inactivate(p1);
			factory->inactivate(p2);
			MovingParticle::setNeighbors(p, pnb[1], pnb[0]);
			MovingParticle::setNeighbors((MovingParticle*)ev.q, qnb[1], qnb[0]);*/
		}
		factory->inactivate(p);
		time = ev.t;
	}
	//printf("counts = %d\n", counts);
	return parts;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "UniformConnectedness", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [X Y] = UniformConnectedness(P, E)");
		return;
	}
	ParticleSimulator simulator;
	//Points
	vector<CParticleF> points; 
	const int* dimsP;
	vector<float> P0;
	mxClassID classIdP;
	int ndimP;
	LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
	for (int i = 0; i < dimsP[0]; ++i)
	{
		float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
		float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
		points.push_back(CParticleF(x, y));
	}
	//edges (indices to the points)
	vector<pair<int, int>> E;
	vector<int> T0;
	mxClassID classIdT;
	int ndimT;
	const int* dimsT;
	LoadData(T0, prhs[1], classIdT, ndimT, &dimsT);
	E = indices2pairs(T0, dimsT);

	simulator.Prepare(points, E);
	vector<Snapshot> parts = Connect();

	if (nlhs >= 1) 
	{
		plhs[0] = Snapshot::StoreSnapshots(parts);
	}
	/*if (nlhs >= 2)
	{
		plhs[1] = simulator.SaveParticles();
	}
	if (nlhs >= 3)
	{
		plhs[2] = Snapshot::StoreSnapshots(simulator.closedRegions);
	}
	if (nlhs >= 4)
	{
		plhs[3] = simulator.SaveDoneEvents();
	}*/

	ParticleFactory::getInstance()->clean();
	mexUnlock();
}

