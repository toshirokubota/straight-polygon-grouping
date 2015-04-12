#include <ParticleSimulator.h>
#include <map>
#include <szMiscOperations.h>
#include <MiscGeometry.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <szmexutilitytemplate.h>

bool
ParticleSimulator::Simulate(float endtime, float delta, bool bdebug)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	float snapTime = delta;
	bool bSuccess = true;
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		(*it)->updateEvent();
	}
	while (time < endtime)
	{
		vector<Snapshot> shots = Snapshot::TakeSnapshot(time); //temporary
		MovingParticle* p = MovingParticle::getNextEvent();
		if (p == NULL) break;
		if (p->getEvent().t > endtime) break;
		if (p->id == 372 || p->id==330)
			p->id += 0;

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
					/*if (closedRegions.size()==1)
					{
					printf("\nSuspicious trace: %d\n", closedRegions.size());
					for (int k = 0; k < tr.size(); ++k)
					{
					printf("%f %f %f %f %d\n",
					tr[k]->p.m_X, tr[k]->p.m_Y, tr[k]->p0.m_X, tr[k]->p0.m_Y, tr[k]->id);
					}
					printf("\nSuspicious closed region: %d\n", closedRegions.size());
					for (int k = 0; k < regions[i].size(); ++k)
					{
					printf("%f %f %f %f %d\n",
					regions[i][k]->p.m_X, regions[i][k]->p.m_Y,
					regions[i][k]->p0.m_X, regions[i][k]->p0.m_Y, regions[i][k]->id);
					}
					printf("\n");
					}*/
				}
				/*Snapshot shot(0.0f, regions[i]);
				if (find(closedRegions.begin(), closedRegions.end(), shot) == closedRegions.end())
				{
				closedRegions.push_back(shot);
				}*/
			}
		}

		doneEvents.push_back(p->getEvent());
		//for (set<MovingParticle*>::iterator it = factory->updateQueue.begin(); it != factory->updateQueue.end(); ++it)
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			(*it)->updateEvent();
		}
		factory->updateQueue.clear();
		if (p->id == 487)
			break;
	}
	return bSuccess;
}

vector<Edge<CParticleF>*>
_trace(Edge<CParticleF>* edge)
{
	Vertex<CParticleF>*first = edge->u;
	vector<Edge<CParticleF>*> path;

	edge->type = Back; //mark it Back once visited.
	path.push_back(edge);
	while (true)
	{
		//choose the first one in the clockwise order
		Edge<CParticleF>* next = NULL;
		float minAngle = 3 * PI;
		//Walk in clock-wise  order.
		Vertex<CParticleF>* u = edge->u;
		Vertex<CParticleF>* v = edge->v;
		if ((int)v->key.m_X == 217 && (int)v->key.m_Y == 107)
		{
			next = NULL;
		}
		for (int i = 0; i < v->aList.size(); ++i)
		{
			Vertex<CParticleF>* w = v->aList[i]->v;
			float ang = GetVisualAngle2(u->key.m_X, u->key.m_Y, w->key.m_X, w->key.m_Y, v->key.m_X, v->key.m_Y);
			if (ang <= 0) ang += 2 * PI;
			if (ang < minAngle)
			{
				minAngle = ang;
				next = v->aList[i];
			}
		}
		if (next->type == Back) break; //Done.

		edge = next;
		edge->type = Back; //mark it Back once visited.
		path.push_back(edge);
	}

	return path;
}

vector<MovingParticle*>
ParticleSimulator::initializePolygon(vector<Edge<CParticleF>*>& edges)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<MovingParticle*> particles;
	for (int i = 0; i < edges.size(); ++i)
	{
		particles.push_back(factory->makeParticle(edges[i]->u->key, Initial, 0.0f));
		if (edges[i]->u->aList.size() <= 1) //leaf node.
		{
			particles.push_back(factory->makeParticle(edges[i]->u->key, Initial, 0.0f));
		}
	}
	//set the neighborhood
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
		MovingParticle::setNeighbors(p, r, q);
	}
	//calculate velocity if necessary
	for (int i = 0; i < particles.size(); ++i)
	{
		MovingParticle* p = particles[i];
		if (p->isInitialized()) continue;
		CParticleF o = p->getP();
		MovingParticle* q = particles[i == particles.size() - 1 ? 0 : i + 1];
		MovingParticle* r = particles[i == 0 ? particles.size() - 1 : i - 1];
		CParticleF b = q->getP();
		CParticleF a = r->getP();
		if (o == a) //leaf node
		{
			float ang = GetVisualDirection(o.m_X, o.m_Y, b.m_X, b.m_Y) + PI / 2.0;
			a.m_X = o.m_X + cos(ang);
			a.m_Y = o.m_Y + sin(ang);
		}
		else if (o == b) //leaf node
		{
			float ang = GetVisualDirection(o.m_X, o.m_Y, a.m_X, a.m_Y) - PI / 2.0;
			b.m_X = o.m_X + cos(ang);
			b.m_Y = o.m_Y + sin(ang);
		}
		float vx, vy;
		calculateBisectorVelocity(a, o, b, vx, vy);
		p->setVelocity(vx, vy);
	}

	return particles;
}

/*
Trace a forrest and create a polygon for each tree.
*/
bool
traceForrest(vector<Vertex<CParticleF>*>& forrest, vector<Edge<CParticleF>*>& edges)
{
	for (int i = 0; i < edges.size(); ++i)
	{
		edges[i]->type = Tr;
	}

	while (true)
	{

		bool bDone = true;
		for (int j = 0; j < edges.size(); ++j)
		{
			if (edges[j]->type == Tr)
			{
				vector<Edge<CParticleF>*> path = _trace(edges[j]);
				/*if (path.size() >= 49)
				{
					for (int k = 0; k < path.size(); ++k)
					{
						printf("%d %3.3f %3.3f %3.3f %3.3f\n", 
							k, path[k]->u->key.m_X, path[k]->u->key.m_Y,
							path[k]->v->key.m_X, path[k]->v->key.m_Y);
					}
				}*/
				vector<MovingParticle*> points = ParticleSimulator::initializePolygon(path);
				bDone = false;
				break;
			}
		}
		if (bDone)
		{
			break;
		}
	}
	return true;
}

bool
ParticleSimulator::Prepare(vector<CParticleF>& points, vector<pair<int, int>>& E, float delta0)
{
	GraphFactory<CParticleF>* factory = GraphFactory<CParticleF>::GetInstance();
	ParticleFactory* pfactory = ParticleFactory::getInstance();

	vector<Vertex<CParticleF>*> vertices;
	for (int i = 0; i < points.size(); ++i)
	{
		vertices.push_back(factory->makeVertex(points[i]));
	}
	vector<Edge<CParticleF>*> edges;
	for (int i = 0; i < E.size(); ++i)
	{
		pair<int, int> idx = E[i];
		Edge<CParticleF>* edge = factory->makeEdge(vertices[idx.first], vertices[idx.second], 1.0f);
		Edge<CParticleF>* edge2 = factory->makeEdge(vertices[idx.second], vertices[idx.first], 1.0f);
		vertices[idx.first]->Add(edge);
		vertices[idx.second]->Add(edge2);
		edges.push_back(edge);
		edges.push_back(edge2);
	}

	traceForrest(vertices, edges);

	time = 0.0f;
	for (set<MovingParticle*>::iterator it = pfactory->activeSet.begin(); it != pfactory->activeSet.end(); ++it)
	{
		(*it)->update(delta0);
	}
	time += delta0;

	/*for (int i = 0; i < pfactory->particles.size(); ++i)
	{
		MovingParticle* p = pfactory->particles[i];
		p->print();
	}*/
	
	return true;
}

#include <mexFileIO.h>
mxArray*
ParticleSimulator::SaveParticles()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	const int dims[] = { factory->particles.size(), ParticleDumpSize };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < dims[0]; ++i)
	{
		MovingParticle* p = factory->particles[i];
		vector<float> v = p->dump2vector();
		for (int j = 0; j < dims[1]; ++j)
		{
			SetData2(F, i, j, dims[0], dims[1], v[j]);
		}
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

mxArray*
ParticleSimulator::SaveDoneEvents()
{
	const int dims[] = { doneEvents.size(), 5 };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < doneEvents.size(); ++i)
	{
		EventStruct ev = doneEvents[i];
		SetData2(F, i, 0, dims[0], dims[1], (float)ev.p->id);
		SetData2(F, i, 1, dims[0], dims[1], (float)ev.q->id);
		SetData2(F, i, 2, dims[0], dims[1], (float)ev.r->id);
		SetData2(F, i, 0, dims[0], dims[1], (float)ev.t);
		SetData2(F, i, 0, dims[0], dims[1], (float)ev.type);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}


/*
a snapshot provides neighborhood information at a particular point in time. We try to reproduce the configuration of the system 
as accurately as possible. We need to consider the following steps.
1. all particles that are created after the time point need to be removed from the system.
2. all particles that are still active at the time point need to be activated.
3. children of a particle need to be erased if they are created after the time point.
*/
bool
ParticleSimulator::Restore(vector<Snapshot>& snapshots)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	float time0 = 0;
	for (int i = 0; i < snapshots.size(); ++i)
	{
		if (snapshots[i].getTime() > time0)
		{
			time0 = snapshots[i].getTime();
		}
	}
	time = time0;
	for (int i = 0; i < factory->particles.size(); ++i)
	{
		MovingParticle* p = factory->particles[i];
		for (int j = 0; j < 2; ++j)
		{
			if (p->children[j] != NULL && p->children[j]->created > time)
			{
				p->children[j] = NULL;
			}
		}
	}

	vector<MovingParticle*> toremove;
	for (int i = 0; i < factory->particles.size(); ++i)
	{
		MovingParticle* p = factory->particles[i];
		if (p->created > time)
		{
			toremove.push_back(p);
		}
	}
	//activate those that are still active.
	vector<MovingParticle*> toactivate;
	for (int i = 0; i < factory->particles.size(); ++i)
	{
		MovingParticle* p = factory->particles[i];
		if (p->bActive == false && p->created < time && p->time > time)
		{
			toactivate.push_back(p);
		}
	}

	//set neighbors
	for (int i = 0; i < snapshots.size(); ++i)
	{
		float t = snapshots[i].getTime();
		for (int j = 0; j < snapshots[i].size(); ++j)
		{
			int id = snapshots[i].get(j);
			int id0 = snapshots[i].get(j == 0 ? snapshots[i].size() - 1 : j - 1);
			int id2 = snapshots[i].get(j == snapshots[i].size() - 1 ? 0 : j + 1);
			MovingParticle* p = factory->get(id);
			if (p == NULL)
			{
				return false;
			}
			p->setNeighbors(p, factory->get(id0), factory->get(id2));
		}
	}

	//delete future particles from memory
	for (int i = 0; i < toremove.size(); ++i)
	{
		factory->remove(toremove[i]);
	}
	//factory->activeSet.clear();
	for (int i = 0; i < toactivate.size(); ++i)
	{
		factory->activeSet.insert(toactivate[i]);
	}
	//recalculate positions
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		MovingParticle* p = *it;
		p->p = p->project(time);
		p->time = time;
	}
	//recalculate events
	for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
	{
		(*it)->event.type == UnknownEvent;
		(*it)->updateEvent();
	}
	return true;
}

bool
ParticleSimulator::LoadParticles(vector<float>& F, const int* dims)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	factory->clean();
	this->time = 0.0f;
	map<int, MovingParticle*> id2particle;
	vector<int> vchild1;
	vector<int> vchild2;
	vector<int> vparent1;
	vector<int> vparent2;
	vector<int> vprev;
	vector<int> vnext;
	vector<int> veq;
	vector<int> ver;
	for (int i = 0; i < dims[0]; ++i)
	{
		int id = (int)GetData2(F, i, 0, dims[0], dims[1], -1.0f);
		float x0 = GetData2(F, i, 1, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float y0 = GetData2(F, i, 2, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float x = GetData2(F, i, 3, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float y = GetData2(F, i, 4, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float vx = GetData2(F, i, 5, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		float vy = GetData2(F, i, 6, dims[0], dims[1], std::numeric_limits<float>::quiet_NaN());
		int pid = (int)GetData2(F, i, 7, dims[0], dims[1], -1.0f);
		int nid = (int)GetData2(F, i, 8, dims[0], dims[1], -1.0f);
		MovingParticleType type = int2ParticleType((int)GetData2(F, i, 9, dims[0], dims[1], -1.0f));
		float created = GetData2(F, i, 10, dims[0], dims[1], 0.0f);
		float time = GetData2(F, i, 11, dims[0], dims[1], 0.0f);
		float ref = GetData2(F, i, 12, dims[0], dims[1], 0.0f);
		bool bActive = GetData2(F, i, 13, dims[0], dims[1], 0.0f)>0.0f ? true : false;
		bool bInitialized = (bool)GetData2(F, i, 14, dims[0], dims[1], 0.0f)>0.0f ? true : false;
		bool bUnstable = (bool)GetData2(F, i, 15, dims[0], dims[1], 0.0f)>0.0f ? true : false;
		int parent1 = (int)GetData2(F, i, 16, dims[0], dims[1], -1.0f);
		int parent2 = (int)GetData2(F, i, 17, dims[0], dims[1], -1.0f);
		int child1 = (int)GetData2(F, i, 18, dims[0], dims[1], -1.0f);
		int child2 = (int)GetData2(F, i, 19, dims[0], dims[1], -1.0f);
		EventType etype = int2EventType((int)GetData2(F, i, 20, dims[0], dims[1], 0.0f));
		int eq = (int)GetData2(F, i, 21, dims[0], dims[1], -1.0f);
		int er = (int)GetData2(F, i, 22, dims[0], dims[1], -1.0f);
		float etime = GetData2(F, i, 23, dims[0], dims[1], 0.0f);
		MovingParticle* p = factory->makeParticle(CParticleF(x0, y0), type, created);
		p->p.m_X = x;
		p->p.m_Y = y;
		p->v[0] = vx;
		p->v[1] = vy;
		p->time = time;
		p->bActive = bActive;
		p->bInitialized = bInitialized;
		p->bUnstable = bUnstable;
		p->event.type = etype;
		p->event.t = etime;
		p->reflexive = ref;
		vparent1.push_back(parent1);
		vparent2.push_back(parent2);
		vchild1.push_back(child1);
		vchild2.push_back(child2);
		p->event.p = p;
		vprev.push_back(pid);
		vnext.push_back(nid);
		veq.push_back(eq);
		ver.push_back(er);
		id2particle[p->id] = p;
		if (p->bActive==false)
		{
			factory->inactivate(p);
		}
		if (p->time > this->time)
		{
			this->time = p->time;
		}
	}
	//now  set prev, next, etc.
	for (int i = 0; i < factory->particles.size(); ++i)
	{
		MovingParticle* p = factory->particles[i];
		p->prev = vprev[i] >= 0 ? id2particle[vprev[i]] : NULL;
		p->next = vnext[i] >= 0 ? id2particle[vnext[i]] : NULL;
		p->event.q = veq[i] >= 0 ? id2particle[veq[i]] : NULL;
		p->event.r = ver[i] >= 0 ? id2particle[ver[i]] : NULL;
		p->parents[0] = vparent1[i] >= 0 ? id2particle[vparent1[i]] : NULL;
		p->parents[1] = vparent2[i] >= 0 ? id2particle[vparent2[i]] : NULL;
		p->children[0] = vchild1[i] >= 0 ? id2particle[vchild1[i]] : NULL;
		p->children[1] = vchild2[i] >= 0 ? id2particle[vchild2[i]] : NULL;
		if (p->id >= MovingParticle::_id)
		{
			MovingParticle::_id = p->id + 1;
		}
	}

	return true;
}

mxArray*
ParticleSimulator::SaveConvexity()
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	const int dims[] = { factory->particles.size(), 3 };
	vector<float> F(dims[0] * dims[1]);
	for (int i = 0; i < dims[0]; ++i)
	{
		MovingParticle* p = factory->particles[i];
		CParticleF pr = p->project(p->created + 0.1);
		SetData2(F, i, 0, dims[0], dims[1], pr.m_X);
		SetData2(F, i, 1, dims[0], dims[1], pr.m_Y);
		SetData2(F, i, 2, dims[0], dims[1], p->reflexive<=0 ? 0.0f: 1.0f);
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

