#include <DistanceMapUtility.h>
#include <szlocalExtrema.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>

vector<CParticleF> 
localMaximaPoints(const vector<float>& dmap, 
			      float threshold,
				  const int* dims)
{
	vector<unsigned char> lmax(dims[0]*dims[1], (unsigned char)0);
	vector<unsigned char> mask(dims[0]*dims[1], (unsigned char)1);
	vector<int> nbh = MakeEightNeighborhood(2, dims);
	//vector<int> nbh = MakeNeighborhood(2, 2, dims);
	LocalMaximum(lmax, dmap, mask, nbh, false, 2, dims);

	vector<CParticleF> locations;
	for(int i=0; i<dims[1]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			if(GetData2(lmax, j, i, dims[0], dims[1], (unsigned char)0))
			{
				float dval = GetData2(dmap, j, i, dims[0], dims[1], (float)0);
				if(dval >threshold)
				{
					locations.push_back(CParticleF(j, i, 0, dval));
				}
			}
		}
	}
	return locations;
}

vector<CParticleFL> 
clusterLocalMaximaPoints(const vector<CParticleF>& locs,
						 const vector<float>& dmap,
						 float ratio, 
						 const int* dims)
{
	vector<Node<CParticleF>*> nodes;
	for(int i=0; i<locs.size(); ++i)
	{
		nodes.push_back(makeset(locs[i]));
	}
	for(int i=0; i<locs.size(); ++i)
	{
		float dval = GetData2(dmap, locs[i].m_X, locs[i].m_Y, dims[0], dims[1], (float)0);
		for(int j=0; j<locs.size(); ++j)
		{
			if(findset(nodes[i]) == findset(nodes[j])) continue;
			if(Distance(locs[i], locs[j]) < ratio * dval)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	
	vector<Node<CParticleF>*> clustered = clusters(nodes);
	vector<CParticleFL> labeled;
	for(int i=0; i<locs.size(); ++i)
	{
		int k = distance(clustered.begin(), find(clustered.begin(), clustered.end(), findset(nodes[i])));
		CParticleFL p(locs[i]);
		p.m_Label = k+1;
		labeled.push_back(p);
	}
	for(int i=0; i<locs.size(); ++i)
	{
		delete nodes[i];
		nodes[i] = NULL;
	}
	return labeled;
}

void
CentricityTransform(vector<float>& C, 
					const vector<float>& D,
					const int* dims)
{
	int xoff[] = {-1, 0, 1, -1, 1, -1, 0, 1};
	int yoff[] = {-1, -1, -1, 0, 0, 1, 1, 1};
	for(int i=0; i<dims[1]; ++i)
	{
		for(int j=0; j<dims[0]; ++j)
		{
			int larger = 0;
			int same = 0;
			float dval = GetData2(D, j, i, dims[0], dims[1], (float)0);
			for(int k=0; k<8; ++k)
			{
				float dval2 = GetData2(D, j+xoff[k], i+yoff[k], dims[0], dims[1], (float)0);
				if(dval2<=0) continue;
				if(dval > dval2) larger++;
				else if(dval==dval2) same++;
			}
			if(same <= 2) 
			{
				larger += same; //permit pixels along a flat ridge
			}
			SetData2(C, j, i, dims[0], dims[1], (float)larger/8);
		}
	}
}

/*
Trace upward along steepest ascent paths.
*/
vector<CParticleFL>
traceUpward(CParticleFL p,
			const vector<float>& dmap,
			const int* dims)
{
	vector<CParticleFL> traced;
	vector<CParticleFL> Q(1, p);
	vector<unsigned char> map(dmap.size(), (unsigned char)0);
	while(Q.empty() == false)
	{
		for(int k=0; k<Q.size(); ++k)
		{
			SetData2(map, Q[k].m_X, Q[k].m_Y, dims[0], dims[1], (unsigned char)1);
		}
		vector<CParticleFL> Q2;
		for(int i=0; i<Q.size(); ++i)
		{
			int x = Q[i].m_X;
			int y = Q[i].m_Y;
			float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
			float dmax = dval;
			vector<CParticleFL> nexts;
			bool bPeak = true;
			for(int k=0; k<NumNeighbors8; ++k)
			{
				int x2 = x + XOffset8[k];
				int y2 = y + YOffset8[k];
				float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
				if(dval2 > dval) bPeak = false;
				if(GetData2(map, x2, y2, dims[0], dims[1], (unsigned char)1)==0)
				{
					if(dval2 >= dmax)
					{
						CParticleFL q(x2, y2, 0, 0, Q[i].m_Label); //propagate the label
						if(dval2 > dmax)
						{
							nexts.clear();
						}
						nexts.push_back(q);
					}
				}
			}
			if(bPeak) Q[i].m_Flag = true;
			traced.push_back(Q[i]);

			for(int k=0; k<nexts.size(); ++k)
			{
				if(find(Q2.begin(), Q2.end(), nexts[k])==Q2.end())
				{
					Q2.push_back(nexts[k]);
				}
			}
		}
		Q = Q2;
	}
	return traced;
}

/*
Trace upward along steepest ascent paths.
Use a representative position for multiple choices.
*/
vector<CParticleF>
traceUpwardRepresentative(CParticleF p,
			const vector<float>& dmap,
			const int* dims)
{
	vector<CParticleF> traced;
	while(true)
	{
		traced.push_back(p);
		int x = p.m_X;
		int y = p.m_Y;
		float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
		float dmax = dval;
		vector<CParticleF> nexts;
		for(int k=0; k<NumNeighbors8; ++k)
		{
			int x2 = x + XOffset8[k];
			int y2 = y + YOffset8[k];
			float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
			if(dval2>dval && dval2 >= dmax)
			{
				CParticleF q(x2, y2, 0, dval2); 
				if(dval2 > dmax)
				{
					nexts.clear();
				}
				nexts.push_back(q);
				dmax = dval2;
			}
		}
		if(nexts.empty()) break;
		else if(nexts.size()==1) p = nexts[0];
		else if(nexts.size()==2)
		{
			if(Abs(nexts[0].m_X-nexts[1].m_X)<=1 && Abs(nexts[0].m_Y-nexts[1].m_Y)<=1)
			{
				p = representativePosition(nexts);
			}
			else
			{
				break; //too ambiguous to committ
			}
		}
		else p = representativePosition(nexts);
	}
	return traced;
}


/*
Trace upward along steepest ascent paths.
In case of tie, increase the neighbors until a single max emerges.
It caps at maxWidth. If reached the max, use a representative.
*/
vector<CParticleF>
traceUpwardVariableNeighbor(CParticleF p,
			const vector<float>& dmap,
			int maxWidth,
			const int* dims)
{
	vector<CParticleF> traced;
	bool bDone = false;
	while(!bDone)
	{
		int x = p.m_X;
		int y = p.m_Y;
		float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
		p.m_Life = dval;
		traced.push_back(p);
		float maxrate = 0;
		float thres = 0.;
		for(int w=1; w<=maxWidth; ++w)
		{
			vector<CParticleF> nexts;
			int countLarger = 0;
			for(int y2=y-w; y2<=y+w; ++y2)
			{
				for(int x2=x-w; x2<=x+w; ++x2)
				{
					if(Abs(x-x2)==w || Abs(y-y2)==w) 
					{
						float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
						if(dval2 > dval) countLarger++;
						float rate = (dval2-dval); // / sqrt((float)(x2-x)*(x2-x) + (y2-y)*(y2-y));
						if(rate>thres && rate >= maxrate)
						{
							CParticleF q(x2, y2, 0, dval2); 
							if(rate > maxrate)
							{
								nexts.clear();
							}
							nexts.push_back(q);
							maxrate = rate;
						}
					}
				}
			}
			if(traced.size()>5 && countLarger < 3) 
			{
				bDone = true;
				break; //try to stop at a ridge.
			}
			if(nexts.empty()) 
			{
				bDone = true;
				break;
			}
			else if(nexts.size()==1)
			{
				p = nexts[0];
				break;
			}
			else if(w == maxWidth)
			{
				p = representativePosition(nexts);
			}
		}
	}
	return traced;
}


/*
Trace upward along steepest ascent paths.
In case of tie, increase the neighbors until a single max emerges.
It caps at maxWidth. If reached the max, use a representative.
*/
vector<CParticleF>
traceUpwardPerptoShape(CParticleF p,
			const vector<float>& dmap,
			const CParticleF& adj,
			const int* dims)
{
	CParticleF shape_dir(p.m_X-adj.m_X, p.m_Y-adj.m_Y);
	shape_dir = UnitVector(shape_dir);

	vector<CParticleF> traced;
	while(true)
	{
		int x = p.m_X;
		int y = p.m_Y;
		float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
		p.m_Life = dval;
		traced.push_back(p);
		float maxrate = 0;
		float thres = 0.0;
		vector<CParticleF> nexts;
		int countLarger = 0;
		for(int k=0; k<NumNeighbors8; ++k)
		{
			int x2 = x + XOffset8[k];
			int y2 = y + YOffset8[k];
			float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
			float rate = (dval2-dval);// / sqrt((float)(x2-x)*(x2-x) + (y2-y)*(y2-y));
			if(dval2 > dval) countLarger++;
			if(rate>thres && rate >= maxrate)
			{
				CParticleF q(x2, y2); 
				if(rate > maxrate)
				{
					nexts.clear();
				}
				nexts.push_back(q);
				maxrate = rate;
			}
		}
		if(traced.size()>5 && countLarger < 3) break; //try to stop at a ridge.

		if(nexts.empty()) break;
		else if(nexts.size()==1)
		{
			p = nexts[0];
		}
		else
		{
			//find a point that gives a trace perpendicular to the shape
			float min_dpval = 1.0;
			for(int i=0; i<nexts.size(); ++i)
			{
				float dx = nexts[i].m_X - x;
				float dy = nexts[i].m_Y - y;
				float dpval = Abs((dx*shape_dir.m_X + dy*shape_dir.m_Y)/sqrt(dx*dx+dy*dy));
				if(dpval < min_dpval)
				{
					min_dpval = dpval;
					p = nexts[i];
				}
			}
		}
	}
	return traced;
}

/*
Trace upward along steepest ascent paths.
In case of tie, it tries to continue the direction it came from.
*/
vector<CParticleF>
traceUpwardPreferStraight(CParticleF p,
						  const vector<float>& dmap,
						  CParticleF prev,
						  const int* dims)
{
	vector<CParticleF> traced;
	while(true)
	{
		int x = p.m_X;
		int y = p.m_Y;
		float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
		p.m_Life = dval;
		traced.push_back(p);
		float maxrate = 0;
		float thres = 0.0;
		vector<CParticleF> nexts;
		int countLarger = 0;
		for(int k=0; k<NumNeighbors8; ++k)
		{
			int x2 = x + XOffset8[k];
			int y2 = y + YOffset8[k];
			float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
			float rate = (dval2-dval);// / sqrt((float)(x2-x)*(x2-x) + (y2-y)*(y2-y));
			if(dval2 > dval) countLarger++;
			if(rate>thres && rate >= maxrate)
			{
				CParticleF q(x2, y2); 
				if(rate > maxrate)
				{
					nexts.clear();
				}
				nexts.push_back(q);
				maxrate = rate;
			}
		}
		if(traced.size()>5 && countLarger < 3) break; //try to stop at a ridge.

		CParticleF next;
		if(nexts.empty()) break;
		else if(nexts.size()==1)
		{
			next = nexts[0];
		}
		else
		{
			//find a point that gives a trace along the previous point
			float dx = x - prev.m_X;
			float dy = y - prev.m_Y;
			float max_dpval = -1.0;
			for(int i=0; i<nexts.size(); ++i)
			{
				float dx2 = nexts[i].m_X - x;
				float dy2 = nexts[i].m_Y - y;
				float dpval = Abs(dx*dx2 + dy*dy2) /(sqrt(dx2*dx2+dy2*dy2) * sqrt(dx*dx + dy*dy));
				if(dpval < max_dpval)
				{
					max_dpval = dpval;
					next = nexts[i];
				}
			}
		}
		prev = p;
		p = next;
	}
	return traced;
}

/*
Trace downward. All paths, not necessarily the steepest one, are considered.
*/
vector<CParticleFL>
traceDownward(const vector<CParticleFL>& seeds,
			  const vector<float>& dmap,
			  const int* dims)
{
	vector<CParticleFL> traced;
	vector<CParticleFL> Q = seeds;
	vector<unsigned char> map(dmap.size(), (unsigned char)0);
	while(Q.empty() == false)
	{
		for(int k=0; k<Q.size(); ++k)
		{
			SetData2(map, Q[k].m_X, Q[k].m_Y, dims[0], dims[1], (unsigned char)1);
		}
		vector<CParticleFL> Q2;
		for(int i=0; i<Q.size(); ++i)
		{
			int x = Q[i].m_X;
			int y = Q[i].m_Y;
			float dval = GetData2(dmap, x, y, dims[0], dims[1], (float)0);
			float dmax = dval;
			vector<CParticleFL> nexts;
			bool bBottom = true;
			for(int k=0; k<NumNeighbors8; ++k)
			{
				int x2 = x + XOffset8[k];
				int y2 = y + YOffset8[k];
				float dval2 = GetData2(dmap, x2, y2, dims[0], dims[1], (float)dval);
				if(dval2 < dval) bBottom = false;
				if(GetData2(map, x2, y2, dims[0], dims[1], (unsigned char)1)==0)
				{
					if(dval2 < dval)
					{
						CParticleFL q(x2, y2, 0, 0, Q[i].m_Label); //propagate the label
						nexts.push_back(q);
					}
				}
			}
			Q[i].m_Flag = bBottom;
			traced.push_back(Q[i]);

			for(int k=0; k<nexts.size(); ++k)
			{
				if(find(Q2.begin(), Q2.end(), nexts[k])==Q2.end())
				{
					Q2.push_back(nexts[k]);
				}
			}
		}
		Q = Q2;
	}
	return traced;
}

/*
a single point seed version.
*/
vector<CParticleFL>
traceDownward(CParticleFL seed,
			  const vector<float>& dmap,
			  const int* dims)
{
	vector<CParticleFL> seeds(1, seed);
	return traceDownward(seeds, dmap, dims);
}
