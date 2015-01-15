#include <Snapshot.h>
#include <szmexutilitytemplate.h>
#include <mexFileIO.h>

vector<Snapshot>
Snapshot::TakeSnapshot(float time)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<Snapshot> polygons;
	set<MovingParticle*> pset;
	while (true)
	{
		bool bdone = true;
		for (set<MovingParticle*>::iterator it = factory->activeSet.begin(); it != factory->activeSet.end(); ++it)
		{
			MovingParticle* p = *it;
			if (pset.find(p) == pset.end())
			{
				vector<MovingParticle*> vp = MovingParticle::vectorize(p);
				Snapshot shot(time);
				for (int i = 0; i < vp.size(); ++i)
				{
					pset.insert(vp[i]);
					shot.add(vp[i]->getId());
				}
				bdone = false;
				polygons.push_back(shot);
			}
		}
		if (bdone) break;
	}
	return polygons;
}

Snapshot
Snapshot::LoadSnapshot(const mxArray* ptr)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	vector<float> F;
	mxClassID classId;
	int ndim;
	const int* dims;
	LoadData(F, ptr, classId, ndim, &dims);
	if (F.size()<4)
	{
		return Snapshot();
	}
	else
	{
		Snapshot shot(-1.0f);
		for (int i = 0; i < dims[0]; ++i)
		{
			float x = GetData2(F, i, 0, dims[0], dims[1], 0.0f);
			float y = GetData2(F, i, 1, dims[0], dims[1], 0.0f);
			float t = GetData2(F, i, 2, dims[0], dims[1], 0.0f);
			int id = (int)GetData2(F, i, 3, dims[0], dims[1], -1.0f);
			if (shot.time < 0)
			{
				shot.time = t;
			}
			else if (t != shot.time)
			{
				mexErrMsgTxt("LoadSnapshot: Inconsistent time stamps are found.");
			}
			shot.add(id);
		}
		return shot;
	}
}

mxArray*
Snapshot::StoreSnapshot(Snapshot& snapshot)
{
	ParticleFactory* factory = ParticleFactory::getInstance();
	const int dims[] = { snapshot.particles.size(), 4 };
	vector<float> F(dims[0] * dims[1]);
	for (int j = 0; j < dims[0]; ++j)
	{
		MovingParticle* p = factory->get(snapshot.particles[j]);
		CParticleF p0 = p->project(Max(p->getCreatedTime(),snapshot.time));
		SetData2(F, j, 0, dims[0], dims[1], p0.m_X);
		SetData2(F, j, 1, dims[0], dims[1], p0.m_Y);
		SetData2(F, j, 2, dims[0], dims[1], snapshot.time);
		SetData2(F, j, 3, dims[0], dims[1], (float)p->getId());
	}
	return StoreData(F, mxSINGLE_CLASS, 2, dims);
}

vector<Snapshot> 
Snapshot::LoadSnapshots(const mxArray* ptr)
{
	int n = mxGetNumberOfElements(ptr);
	vector<Snapshot> snapshots(n);
	for (int i = 0; i < n; ++i)
	{
		mxArray* cell = mxGetCell(ptr, i);
		snapshots[i] = LoadSnapshot(cell);
	}
	return snapshots;
}


mxArray*
Snapshot::StoreSnapshots(vector<Snapshot>& snapshots)
{
	mxArray* plhs;
	const int dims[] = { snapshots.size(), 1 };
	plhs = mxCreateCellArray((mwSize)2, (const mwSize*)dims);
	for (int i = 0; i<snapshots.size(); ++i)
	{
		mxArray* a = StoreSnapshot(snapshots[i]);
		mxSetCell(plhs, i, a);
	}
	return plhs;
}