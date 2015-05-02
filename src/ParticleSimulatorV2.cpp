#include <ParticleSimulatorV2.h>
#include <map>
#include <szMiscOperations.h>
#include <MiscGeometry.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <szmexutilitytemplate.h>
#include <unordered_set>

bool
ParticleSimulatorV2::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		(*it)->updateEvent();
	}
	int iter = 0;
	while (time < endtime)
	{
		iter++;
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time); //temporary
		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;
		if (p->getEvent().t > endtime) break;

		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			(*it)->update(p->getEvent().t - time);
		}
		time = p->getEvent().t;
		if (p->applyEvent() == false) break;
		p->getEvent().print();

		MovingParticle::removeUnstable();
		if (time >= snapTime)
		{
			vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			snapTime += delta;
		}
		MovingParticle::quickFinish();
		if (bdebug && MovingParticle::sanityCheck() == false)
		{
			printf("Violation of sanity check found at %f.\n", time);
			vector<Snapshot> shots = Snapshot::TakeSnapshot(time);
			snapshots.insert(snapshots.end(), shots.begin(), shots.end());
			bSuccess = false;
			break;
		}

		//find closed regions
		vector<vector<MovingParticle*>> regions = MovingParticle::clusterParticles();
		unordered_set<CParticleF,hash<CParticleF>> loopy;
		for (int i = 0; i < regions.size(); ++i)
		{
			//if (clockWise(regions[i]) > 0)
			{
				vector<MovingParticle*> tr = MovingParticle::traceBackPolygon(regions[i]);
				vector<vector<MovingParticle*>> areas = MovingParticle::closedRegions(tr);
				for (int j = 0; j < areas.size(); ++j)
				{
					Snapshot shot(0.0f, areas[j]);
					if (find(closedRegions.begin(), closedRegions.end(), shot) == closedRegions.end())
					{
						closedRegions.push_back(shot);
					}
					for (int k = 0; k < areas[j].size(); ++k) {
						loopy.insert(areas[j][k]->getP0());
					}
				}
			}
		}

		//find descendents of loopy particles and inactivate them.
		vector<MovingParticle*> vp2inactivate;
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			vector<MovingParticle*> vp(1, *it);
			vector<MovingParticle*> tr = MovingParticle::traceBackPolygon(vp);
			for (int k = 0; k < tr.size(); ++k) 
			{
				if (loopy.find(tr[k]->getP0()) != loopy.end()) 
				{
					vp2inactivate.push_back(vp[0]);
					break;
				}
			}
		}
		for (int i = 0; i < vp2inactivate.size(); ++i)
		{
			factory->inactivate(vp2inactivate[i]);
		}

		doneEvents.push_back(p->getEvent());
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			(*it)->updateEvent();
		}
		factory->updateQueue.clear();
	}
	return bSuccess;
}

