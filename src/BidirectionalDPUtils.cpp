#include <BidirectionalDPUtils.h>
#include <szIndexedData.h>
#include <szMiscOperations.h>
#include <szMexUtility.h>
#include <DisjointSet.h>

void
PrintCycle(vector<Vertex<FragmentInfo>*>& cycle)
{
	if(cycle.empty()) 
	{
		printf("PrintCycle: empty cycle.\n");
		return;
	}

	for(int i=0; i<cycle.size(); ++i)
	{
		printf("%d: [%d,%d] a=%3.3f\td=%3.3f\tf=%3.3f\t%d\n", i+1, cycle[i]->key.contourID+1, cycle[i]->key.pointID+1,
			cycle[i]->key.angle, cycle[i]->d, cycle[i]->f, cycle[i]->color);
	}
	printf("Path cost = %f\n", cycle[cycle.size()-1]->d);
}


/*
Stitch a number of EQW contours together into one EQW contour.
The orientation of input contours are assumed to be arranged properly in
by splitContours().
*/
ContourEQW
stitchContour(vector<Vertex<FragmentInfo>*>& cycle, 
			  vector<ContourEQW>& contours)
{
	ContourEQW c;
	for(int i=0; i<cycle.size(); ++i)
	{
		int id1 = cycle[i]->key.contourID;
		int p1 = cycle[i]->key.pointID;
		int id2 = cycle[(i+1)%cycle.size()]->key.contourID;
		int p2 = cycle[(i+1)%cycle.size()]->key.pointID;
		if(id1==id2)
		{
			if(p1<p2)
			{
				c.A.insert(c.A.end(), contours[id1].A.begin()+p1, contours[id1].A.begin()+p2-1);
				c.B.insert(c.B.end(), contours[id1].B.begin()+p1, contours[id1].B.begin()+p2-1);
				c.C.insert(c.C.end(), contours[id1].C.begin()+p1, contours[id1].C.begin()+p2-1);
				c.D.insert(c.D.end(), contours[id1].D.begin()+p1, contours[id1].D.begin()+p2-1);
				c.X.insert(c.X.end(), contours[id1].X.begin()+p1, contours[id1].X.begin()+p2-1);
				c.Y.insert(c.Y.end(), contours[id1].Y.begin()+p1, contours[id1].Y.begin()+p2-1);
				c.Strength.insert(c.Strength.end(), contours[id1].Strength.begin()+p1, contours[id1].Strength.begin()+p2);
			}
			else if(p1>p2)
			{
				for(int j=p1; j>p2; j--)
				{
					c.A.push_back(contours[id1].A[j]);
					c.B.push_back(contours[id1].B[j]);
					c.C.push_back(contours[id1].C[j]);
					c.D.push_back(contours[id1].D[j]);
					c.X.push_back(contours[id1].X[j]);
					c.Y.push_back(contours[id1].Y[j]);
					c.Strength.push_back(contours[id1].Strength[j]);
				}
			}
			else //must be source and destination
			{
				c.A.push_back(contours[id1].A[p1]);
				c.B.push_back(contours[id1].B[p1]);
				c.C.push_back(contours[id1].C[p1]);
				c.D.push_back(contours[id1].D[p1]);
				c.X.push_back(contours[id1].X[p1]);
				c.Y.push_back(contours[id1].Y[p1]);
				c.Strength.push_back(contours[id1].Strength[p1]);
			}
		}
		else
		{
			c.A.push_back(contours[id1].A[p1]);
			c.B.push_back(contours[id1].B[p1]);
			c.C.push_back(contours[id1].C[p1]);
			c.D.push_back(contours[id1].D[p1]);
			c.X.push_back(contours[id1].X[p1]);
			c.Y.push_back(contours[id1].Y[p1]);
			c.Strength.push_back(contours[id1].Strength[p1]);
		}
	}
	return c;
}

vector<ContourEQW> 
splitContoursBySaliency(vector<ContourEQW>& contours, 
						float thres)
{
	vector<ContourEQW> contours2;
	for(int i=0; i<contours.size(); ++i)
	{
		int j=1;
		while(j<contours[i].size()-1)
		{
			if(contours[i].Strength[j-1] < thres && contours[i].Strength[j+1] < thres)
			{
				if(contours[i].Strength[j] > thres)
				{
					contours2.push_back(contours[i].extract(0, j-1));
				}
				else
				{
					contours2.push_back(contours[i].extract(0, j));
				}
			}
			else if(contours[i].Strength[j-1] > thres && contours[i].Strength[j+1] < thres)
			{
				if(contours[i].Strength[j] > thres)
				{
					contours2.push_back(contours[i].extract(0, j));
				}
				else
				{
					contours2.push_back(contours[i].extract(0, j-1));
				}
			}
			j++;
		}
	}
	return contours2;
}

/*
In addition to the length requirement, this version also splits a contour at the point closest to the start position.
*/
vector<ContourEQW> 
splitContoursByLengthAndAtSource(vector<ContourEQW>& contours0, 
								int maxLen, 
								CParticleF& start)
{
	vector<ContourEQW> contours = splitContours(contours0, maxLen);
	float mind = std::numeric_limits<float>::infinity();
	FragmentInfo info;
	for(int i=0; i<contours.size(); ++i)
	{
		for(int j=0; j<contours[i].size(); ++j)
		{
			float d = Distance(contours[i].X[j], contours[i].Y[j], start.m_X, start.m_Y);
			if(d < mind)
			{
				mind = d;
				info.contourID = i;
				info.pointID = j;
			}
		}
	}
	assert(mind < std::numeric_limits<float>::infinity());
	if(info.pointID>0 && info.pointID < contours[info.contourID].size()-1)
	{
		ContourEQW cnew;
		cnew.A.insert(cnew.A.begin(), contours[info.contourID].A.begin()+info.pointID, contours[info.contourID].A.end());
		cnew.B.insert(cnew.B.begin(), contours[info.contourID].B.begin()+info.pointID, contours[info.contourID].B.end());
		cnew.C.insert(cnew.C.begin(), contours[info.contourID].C.begin()+info.pointID, contours[info.contourID].C.end());
		cnew.D.insert(cnew.D.begin(), contours[info.contourID].D.begin()+info.pointID, contours[info.contourID].D.end());
		cnew.X.insert(cnew.X.begin(), contours[info.contourID].X.begin()+info.pointID, contours[info.contourID].X.end());
		cnew.Y.insert(cnew.Y.begin(), contours[info.contourID].Y.begin()+info.pointID, contours[info.contourID].Y.end());
		cnew.Strength.insert(cnew.Strength.begin(), contours[info.contourID].Strength.begin()+info.pointID, contours[info.contourID].Strength.end());

		contours[info.contourID].A.erase(contours[info.contourID].A.begin()+info.pointID, contours[info.contourID].A.end());
		contours[info.contourID].B.erase(contours[info.contourID].B.begin()+info.pointID, contours[info.contourID].B.end());
		contours[info.contourID].C.erase(contours[info.contourID].C.begin()+info.pointID, contours[info.contourID].C.end());
		contours[info.contourID].D.erase(contours[info.contourID].D.begin()+info.pointID, contours[info.contourID].D.end());
		contours[info.contourID].X.erase(contours[info.contourID].X.begin()+info.pointID, contours[info.contourID].X.end());
		contours[info.contourID].Y.erase(contours[info.contourID].Y.begin()+info.pointID, contours[info.contourID].Y.end());
		contours[info.contourID].Strength.erase(contours[info.contourID].Strength.begin()+info.pointID, contours[info.contourID].Strength.end());

		contours.insert(contours.begin()+info.contourID, cnew);
	}
	return contours;
}


/*
Collect a number of contours according to the order given by the cycle.
It keeps multiple contours, not like stitchContour in which the contours
are combined into one.
*/
vector<ContourEQW>
collectContours(vector<Vertex<FragmentInfo>*>& cycle, 
			  vector<ContourEQW>& contours)
{
	vector<ContourEQW> vc;
	for(int i=0; i<cycle.size(); ++i)
	{
		ContourEQW c;
		int id1 = cycle[i]->key.contourID;
		int p1 = cycle[i]->key.pointID;
		int id2 = cycle[(i+1)%cycle.size()]->key.contourID;
		int p2 = cycle[(i+1)%cycle.size()]->key.pointID;
		if(id1==id2)
		{
			if(p1<p2)
			{
				c.A.insert(c.A.end(), contours[id1].A.begin()+p1, contours[id1].A.begin()+p2-1);
				c.B.insert(c.B.end(), contours[id1].B.begin()+p1, contours[id1].B.begin()+p2-1);
				c.C.insert(c.C.end(), contours[id1].C.begin()+p1, contours[id1].C.begin()+p2-1);
				c.D.insert(c.D.end(), contours[id1].D.begin()+p1, contours[id1].D.begin()+p2-1);
				c.X.insert(c.X.end(), contours[id1].X.begin()+p1, contours[id1].X.begin()+p2-1);
				c.Y.insert(c.Y.end(), contours[id1].Y.begin()+p1, contours[id1].Y.begin()+p2-1);
				c.Strength.insert(c.Strength.end(), contours[id1].Strength.begin()+p1, contours[id1].Strength.begin()+p2);
			}
			else if(p1>p2)
			{
				for(int j=p1; j>p2; j--)
				{
					c.A.push_back(contours[id1].A[j]);
					c.B.push_back(contours[id1].B[j]);
					c.C.push_back(contours[id1].C[j]);
					c.D.push_back(contours[id1].D[j]);
					c.X.push_back(contours[id1].X[j]);
					c.Y.push_back(contours[id1].Y[j]);
					c.Strength.push_back(contours[id1].Strength[j]);
				}
			}
			else //must be source and destination
			{
				c.A.push_back(contours[id1].A[p1]);
				c.B.push_back(contours[id1].B[p1]);
				c.C.push_back(contours[id1].C[p1]);
				c.D.push_back(contours[id1].D[p1]);
				c.X.push_back(contours[id1].X[p1]);
				c.Y.push_back(contours[id1].Y[p1]);
				c.Strength.push_back(contours[id1].Strength[p1]);
			}
			i++; //additional increment
		}
		else
		{ //a contour with a single point
			c.A.push_back(contours[id1].A[p1]);
			c.B.push_back(contours[id1].B[p1]);
			c.C.push_back(contours[id1].C[p1]);
			c.D.push_back(contours[id1].D[p1]);
			c.X.push_back(contours[id1].X[p1]);
			c.Y.push_back(contours[id1].Y[p1]);
			c.Strength.push_back(contours[id1].Strength[p1]);
		}
		vc.push_back(c);
	}
	return vc;
}

/*
Split contours so that each length is at most the second argument (length).
*/
vector<ContourEQW> splitContours(vector<ContourEQW>& contours, 
								 int length)
{
	vector<ContourEQW> splitted;
	for(int i=0; i<contours.size(); ++i)
	{
		int m = contours[i].size()/length;
		int rem = contours[i].size() % length;
		m++;
		//figure out how to split a contour of length m 
		vector<int> vlength(m, 0);
		int dif = contours[i].size() - m * length;
		float inc = (float)dif / (float)vlength.size();
		float sum = 0;
		int total = 0;
		for(int k=0; k<vlength.size(); ++k)
		{
			sum += inc;
			vlength[k] = length + (int)sum;
			sum -= (int)sum;
			total += vlength[k];
		}
		if(total > contours[i].size()) //this case may happen due to rounding error
		{
			vlength[0] -= (total - contours[i].size());
		}

		int acc = 0;
		float gap = 1.0e-3;
		for(int j=0; j<m; ++j)
		{
			ContourEQW c = extractContourFragment(contours[i], acc, acc+vlength[j]-1);
			if(c.size()>=2)
			{
				//add a little gap to the end points so that we will have different view angle at the shared point.
				int sz = c.size();
				c.X[0] += (c.X[1]>c.X[0]) ? gap: -gap;
				c.Y[0] += (c.Y[1]>c.Y[0]) ? gap: -gap;
				c.X[sz-1] += (c.X[sz-2]>c.X[sz-1]) ? gap: -gap;
				c.Y[sz-1] += (c.Y[sz-2]>c.Y[sz-1]) ? gap: -gap;
			}
			splitted.push_back(c);
			acc += vlength[j];
		}
	}
	return splitted;
}

/*
Convert a path given by the first argument (vertices) and the second (contours) into a matrix and
return the result in the third argument.
It also returns the number of columns as the result in the fourth argument (ncols).
*/
void
Path2Mat(vector<Vertex<FragmentInfo>*>& vertices,
		  vector<ContourEQW>& contours,
		  vector<float>& res, 
		  int& ncols)
{
	vector<float> mat;
	int n = vertices.size();
	ncols = 14;
	for(int i=0; i<n; ++i)
	{
		Vertex<FragmentInfo>* u = vertices[i];
		Vertex<FragmentInfo>* v = u->pi;
		if(v)
		{
			float x1 = contours[u->key.contourID].X[u->key.pointID];
			float y1 = contours[u->key.contourID].Y[u->key.pointID];
			float x2 = contours[v->key.contourID].X[v->key.pointID];
			float y2 = contours[v->key.contourID].Y[v->key.pointID];

			mat.push_back(x1);
			mat.push_back(y1);
			mat.push_back(u->d);
			mat.push_back(u->f);
			mat.push_back(u->key.angle); //one based index
			mat.push_back(u->key.contourID+1); //one based index
			mat.push_back(u->key.pointID+1); //one based index
			mat.push_back(x2);
			mat.push_back(y2);
			mat.push_back(v->d);
			mat.push_back(v->f);
			mat.push_back(v->key.angle); //one based index
			mat.push_back(v->key.contourID+1); //one based index
			mat.push_back(v->key.pointID+1); //one based index
		}
	}
	res.resize(mat.size()); //transpose the result
	int nrows = mat.size()/ncols;
	for(int i=0; i<nrows; ++i)
	{
		for(int j=0; j<ncols; ++j)
		{
			res[j*nrows+i] = mat[i*ncols+j];
		}
	}
}

/*
THIS FUNCTION IS MOVED TO BidirectionalDPUtilsGraph.cpp.

Generate a set of vertices from a set of contours.
The number of vertices is twice the number of contours (one for the first and the other for the last
point in a contour).
*/
/*vector<Vertex<FragmentInfo>*>
generateVertices(vector<ContourEQW>& contours)
{
	GraphFactory<FragmentInfo>* factory = GraphFactory<FragmentInfo>::GetInstance();
	vector<Vertex<FragmentInfo>*> vertices;
	for(int i=0; i<contours.size(); ++i)
	{
		Vertex<FragmentInfo>* u = factory->makeVertex(FragmentInfo(i, 0));
		Vertex<FragmentInfo>* v = factory->makeVertex(FragmentInfo(i, contours[i].size()-1));
		vertices.push_back(u);
		vertices.push_back(v);
	}
	return vertices;
}*/

/*
Order vertices (1st arg) based on the visual angle seen from a click point (3rd arg).
The 2nd argument carries contour information.
*/
vector<Vertex<FragmentInfo>*>
orderVertices(vector<Vertex<FragmentInfo>*>& vertices,
			  vector<ContourEQW>& contours, 
			  CParticleF& click)
{
	vector<Vertex<FragmentInfo>*> ordered;
	vector<indexedData> vdata;
	for(int i=0; i<vertices.size(); ++i)
	{
		int c = vertices[i]->key.contourID;
		int p = vertices[i]->key.pointID;
		float x = contours[c].X[p];
		float y = contours[c].Y[p];
		indexedData d;
		d.data = GetVisualDirection(x, y, click.m_X, click.m_Y);
		if(d.data < 0)
		{
			 d.data += 2*PI;
		}
		d.index = i;
		vdata.push_back(d);
		vertices[i]->key.angle = d.data;
	}
	sort(vdata.begin(), vdata.end());
	for(int i=0; i<vdata.size(); ++i)
	{
		ordered.push_back(vertices[vdata[i].index]);
	}
	return ordered;
}

/*
Circularly shift vertices so that the source becomes the first one. 
Add a copy of the source at the end as the destination vertex.
*/
vector<Vertex<FragmentInfo>*>
arrangeVertices(vector<Vertex<FragmentInfo>*>& vertices,
				Vertex<FragmentInfo>* source)
{
	int k = distance(vertices.begin(), find(vertices.begin(), vertices.end(), source));
	assert(k < vertices.size());
	vector<Vertex<FragmentInfo>*> result;
	float angle = source->key.angle;
	for(int i=0; i<vertices.size(); ++i)
	{
		int j = (i + k) % vertices.size();
		vertices[j]->key.angle -= angle;
		if(vertices[j]->key.angle < 0)
		{
			vertices[j]->key.angle += 2*PI;
		}
		result.push_back(vertices[j]);
	}
	//GraphFactory<FragmentInfo>* factory = GraphFactory<FragmentInfo>::GetInstance();
	//Vertex<FragmentInfo>* dest = factory->makeVertex(source->key);
	//result.push_back(dest);
	return result;
}

/*
Given a position (x, y) as the 3rd and 4th arguments, find a contour point
that is the closest, and returns the vertex that is closest to the point.
*/
Vertex<FragmentInfo>* 
findSource(vector<Vertex<FragmentInfo>*>& vertices, 
		   vector<ContourEQW>& contours, 
		   float x, float y)
{
	Vertex<FragmentInfo>* source = NULL;
	float d = numeric_limits<float>::infinity();
	FragmentInfo info;
	for(int i=0; i<contours.size(); ++i)
	{
		for(int j=0; j<contours[i].size(); ++j)
		{
			float x2=contours[i].X[j];
			float y2=contours[i].Y[j];
			float d2 = (x-x2)*(x-x2)+(y-y2)*(y-y2);
			if(d2 < d)
			{
				d = d2;
				info.contourID = i;
				info.pointID = j;
			}
		}
	}
	if(d < numeric_limits<float>::infinity())
	{
		int mindif = -1;
		for(int i=0; i<vertices.size(); ++i)
		{
			if(vertices[i]->key.contourID == info.contourID)
			{
				if(mindif<0 || mindif > Abs(info.pointID - vertices[i]->key.pointID))
				{
					mindif = Abs(info.pointID - vertices[i]->key.pointID);
					source = vertices[i];
				}
			}
		}
	}
	/*printf("Source Node = [%d:%d] @ (%2.2f, %2.2f).\n", source->key.contourID+1, source->key.pointID+1,
		contours[source->key.contourID].X[source->key.pointID],
		contours[source->key.contourID].Y[source->key.pointID]);*/
	return source;
}

/*
Perform update on an angularly ordered graph (given as the 1st argument). The 2nd argument
carries the necessary contour information.
The third argument dictates the angular direction (true for contour clock wise, and false for
clock wise).
Currently, the finalized time is set to the iteration number of the update. Thus, the 4th argument
is necessary.
*/
bool 
UpdateDistance(vector<Vertex<FragmentInfo>*>& vertices, 
			   vector<ContourEQW>& contours, 
			   bool orientation,
			   float angleThres,
			   int iter)
{
	bool bchanged = false;
	Vertex<FragmentInfo>* dest = vertices[vertices.size()-1];
	if(orientation)
	{
		for(int i=0; i<vertices.size(); ++i) 
		{
			int cid1 = vertices[i]->key.contourID;
			int pid1 = vertices[i]->key.pointID;
			float d1 = vertices[i]->d;
			for(int j=0; j<vertices.size(); ++j) //include the source
			{
				if(i == j) continue;
				if(vertices[i]->key.angle < vertices[j]->key.angle && Abs(vertices[i]->key.angle-vertices[j]->key.angle)<angleThres)
				{
					int cid2 = vertices[j]->key.contourID;
					int pid2 = vertices[j]->key.pointID;
					float d2 = d1;
					if(cid1 != cid2)
					{
						float d0 = Distance(contours[cid1].X[pid1], contours[cid1].Y[pid1], contours[cid2].X[pid2], contours[cid2].Y[pid2]);
						d2 = d2 + EdgeWeight(d0);
					}
					if(d2 < vertices[j]->d)
					{
						vertices[j]->d = d2;
						vertices[j]->pi = vertices[i];
						vertices[j]->f = iter;
						bchanged = true;
					}
				}
			}
		}
	}
	else
	{
		for(int i=0; i<vertices.size(); ++i) 
		{
			int cid1 = vertices[i]->key.contourID;
			int pid1 = vertices[i]->key.pointID;
			float d1 = vertices[i]->d;
			for(int j=0; j<vertices.size(); ++j) //include the source
			{
				if(i == j) continue;
				if(vertices[i]->key.angle >= vertices[j]->key.angle && Abs(vertices[i]->key.angle-vertices[j]->key.angle)<angleThres)
				{
					int cid2 = vertices[j]->key.contourID;
					int pid2 = vertices[j]->key.pointID;
					float d2 = d1;
					if(cid1 != cid2)
					{
						float d0 = Distance(contours[cid1].X[pid1], contours[cid1].Y[pid1], contours[cid2].X[pid2], contours[cid2].Y[pid2]);
						d2 = d2 + EdgeWeight(d0);
					}
					if(d2 < vertices[j]->d)
					{
						vertices[j]->d = d2;
						vertices[j]->pi = vertices[i];
						vertices[j]->f = iter;
						bchanged = true;
					}
				}
			}
		}
	}
	return bchanged;
}


/*
A slight modification of UpdateDistance.
This version allows non-zero weights on solid fragments. The feature is used to give preference to fragments that are selected by the user.
It is used in BidirectionalDPMPv3.cpp.
*/
bool 
UpdateDistanceWithPreference(vector<Vertex<FragmentInfo>*>& vertices, vector<ContourEQW>& contours, 
							 bool orientation,
							 int iter)
{
	float maxAng = PI/4;
	bool bchanged = false;
	int L = vertices.size()/2;
	Vertex<FragmentInfo>* dest = vertices[vertices.size()-1];
	if(orientation)
	{
		for(int i=0; i<vertices.size(); ++i) 
		{
			int cid1 = vertices[i]->key.contourID;
			int pid1 = vertices[i]->key.pointID;
			for(int j=0; j<i; ++j) //include the source
			{
				if(i-j > L)
				{
					continue;
				}
				int cid2 = vertices[j]->key.contourID;
				int pid2 = vertices[j]->key.pointID;
				float d = vertices[j]->d;
				if(cid1 != cid2)
				{
					float d0 = Distance(contours[cid1].X[pid1], contours[cid1].Y[pid1], contours[cid2].X[pid2], contours[cid2].Y[pid2]);
					d = d + EdgeWeight(d0);
				}
				else //this is where we used fragment specific weight
				{
					float d0 = vertices[i]->key.saliency;
					d += d0;
					/*if(d0<0)
					{
						printf("Negative strength: %f @[%d, %d] - [%d, %d]\n", d0, cid1, pid1, cid2, pid2);
					}*/
				}
				if(d < vertices[i]->d)
				{
					vertices[i]->d = d;
					vertices[i]->pi = vertices[j];
					vertices[i]->f = iter;
					if(cid1 == cid2)
					{
						if(pid1 < pid2)
						{
							for(int k=pid1; k<pid2; ++k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
						else
						{
							for(int k=pid1; k>pid2; --k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
					}
					else
					{
						contours[cid1].Strength[pid1] = (float)iter;
					}
					bchanged = true;
				}
			}
		}
	}
	else
	{
		for(int i=vertices.size()-1; i>=0; --i) //exclude the source
		{
			int cid1 = vertices[i]->key.contourID;
			int pid1 = vertices[i]->key.pointID;
			for(int j=i+1; j<vertices.size(); ++j) //include the dest 
			{
				if(j-i>L)
				{
					continue;
				}
				int cid2 = vertices[j]->key.contourID;
				int pid2 = vertices[j]->key.pointID;
				float d = vertices[j]->d;
				if(cid1 != cid2)
				{
					float d0 = Distance(contours[cid1].X[pid1], contours[cid1].Y[pid1], contours[cid2].X[pid2], contours[cid2].Y[pid2]);
					d = d + EdgeWeight(d0);
				}
				else //this is where we used fragment specific weight
				{
					float d0 = vertices[i]->key.saliency;
					d += d0;
					/*if(d0<0)
					{
						printf("Negative strength (R): %f @[%d, %d] - [%d, %d]\n", d0, cid1, pid1, cid2, pid2);
					}*/
				}
				if(d < vertices[i]->d)
				{
					vertices[i]->d = d;
					vertices[i]->pi = vertices[j];
					vertices[i]->f = iter;
					if(cid1 == cid2)
					{
						if(pid1 < pid2)
						{
							for(int k=pid1; k<pid2; ++k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
						else
						{
							for(int k=pid1; k>pid2; --k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
					}
					else
					{
						contours[cid1].Strength[pid1] = (float)iter;
					}
					bchanged = true;
				}
			}
		}
	}
	return bchanged;
}


/*
This is a slight modification to UpdateDistance. This version can exclude some nodes from being used
in the shortest path.
The fourth argument (excluded) supplies a set of vertices that need to be excluded from the update.
*/
bool 
UpdateDistanceWithExclusions(vector<Vertex<FragmentInfo>*>& vertices, vector<ContourEQW>& contours, 
			   bool orientation,
			   int iter,
			   vector<FragmentInfo>& excluded)
{
	float maxAng = PI/4;
	bool bchanged = false;
	int L = vertices.size()/2;
	Vertex<FragmentInfo>* dest = vertices[vertices.size()-1];
	if(orientation)
	{
		for(int i=0; i<vertices.size(); ++i) 
		{
			int cid1 = vertices[i]->key.contourID;
			int pid1 = vertices[i]->key.pointID;
			for(int j=0; j<i; ++j) //include the source
			{
				if(i-j > L)
				{
					continue;
				}
				if(find(excluded.begin(), excluded.end(), vertices[j]->key) != excluded.end())
				{
					continue;
				}
				int cid2 = vertices[j]->key.contourID;
				int pid2 = vertices[j]->key.pointID;
				float d = vertices[j]->d;
				if(cid1 != cid2)
				{
					float d0 = Distance(contours[cid1].X[pid1], contours[cid1].Y[pid1], contours[cid2].X[pid2], contours[cid2].Y[pid2]);
					d = d + EdgeWeight(d0);
				}
				if(d < vertices[i]->d)
				{
					vertices[i]->d = d;
					vertices[i]->pi = vertices[j];
					vertices[i]->f = iter;
					if(cid1 == cid2)
					{
						if(pid1 < pid2)
						{
							for(int k=pid1; k<pid2; ++k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
						else
						{
							for(int k=pid1; k>pid2; --k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
					}
					else
					{
						contours[cid1].Strength[pid1] = (float)iter;
					}
					bchanged = true;
				}
			}
		}
	}
	else
	{
		for(int i=vertices.size()-1; i>=0; --i) //exclude the source
		{
			int cid1 = vertices[i]->key.contourID;
			int pid1 = vertices[i]->key.pointID;
			for(int j=i+1; j<vertices.size(); ++j) //include the dest 
			{
				if(j-i>L)
				{
					continue;
				}
				if(find(excluded.begin(), excluded.end(), vertices[j]->key) != excluded.end())
				{
					continue;
				}
				int cid2 = vertices[j]->key.contourID;
				int pid2 = vertices[j]->key.pointID;
				float d = vertices[j]->d;
				if(cid1 != cid2)
				{
					float d0 = Distance(contours[cid1].X[pid1], contours[cid1].Y[pid1], contours[cid2].X[pid2], contours[cid2].Y[pid2]);
					d = d + EdgeWeight(d0);
				}
				if(d < vertices[i]->d)
				{
					vertices[i]->d = d;
					vertices[i]->pi = vertices[j];
					vertices[i]->f = iter;
					if(cid1 == cid2)
					{
						if(pid1 < pid2)
						{
							for(int k=pid1; k<pid2; ++k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
						else
						{
							for(int k=pid1; k>pid2; --k)
							{
								contours[cid1].Strength[k] = (float)iter;
							}
						}
					}
					else
					{
						contours[cid1].Strength[pid1] = (float)iter;
					}
					bchanged = true;
				}
			}
		}
	}
	return bchanged;
}

float
getPathCost(vector<Vertex<FragmentInfo>*>& path,
			vector<ContourEQW>& contours)
{
	float cost = 0;
	for(int i=0; i<path.size()-1; ++i)
	{
		int cid1 = path[i]->key.contourID;
		int pid1 = path[i]->key.pointID;
		int cid2 = path[i+1]->key.contourID;
		int pid2 = path[i+1]->key.pointID;
		if(cid1 != cid2)
		{
			float x1 = contours[cid1].X[pid1];
			float y1 = contours[cid1].Y[pid1];
			float x2 = contours[cid2].X[pid2];
			float y2 = contours[cid2].Y[pid2];
			cost = cost + EdgeWeight(Distance(x1, y1, x2, y2));
		}
	}
	return cost;
}

/*
Given an angularly ordered graph, trace a path from the 2nd argument (dest).
*/
vector<Vertex<FragmentInfo>*> 
findPath(vector<Vertex<FragmentInfo>*>& vertices, 
		  Vertex<FragmentInfo>* dest)
{

	vector<Vertex<FragmentInfo>*> cycle;
	Vertex<FragmentInfo>* p = dest;
	while(p)
	{
		//cycle.push_back(p);
		if(find(cycle.begin(), cycle.end(), p)==cycle.end())
		{
			cycle.insert(cycle.begin(), p);
			p = p->pi;
		}
		else
		{
			printf("findPath: a cycle is found.\n");
			for(int i=0; i<cycle.size(); ++i)
			{
				printf("[%d, %d] %f\n", cycle[i]->key.contourID+1, cycle[i]->key.pointID+1, cycle[i]->d);
			}
			break;
		}
	}
	//cycle.erase(cycle.end()-1); //remove the destination.
	for(int i=0; i<cycle.size(); ++i)
	{
		//printf("%d: [%d, %d], %f %d\n", i, cycle[i]->key.contourID+1, cycle[i]->key.pointID+1, cycle[i]->d, (int)cycle[i]->f);
	}
	return cycle;
}

/*
Given an angularly ordered tree (vertices), find leaf vertices.
*/
vector<Vertex<FragmentInfo>*>
findEndPoints(vector<Vertex<FragmentInfo>*>& vertices)
{
	vector<bool> flags(vertices.size(), false);
	for(int i=0; i<vertices.size(); ++i)
	{
		if(vertices[i]->pi)
		{
			int k = distance(vertices.begin(), find(vertices.begin(), vertices.end(), vertices[i]->pi));
			flags[k] = true;
		}
	}
	vector<Vertex<FragmentInfo>*> ends;
	for(int i=0; i<vertices.size(); ++i)
	{
		if(flags[i]==false)
		{
			ends.push_back(vertices[i]);
		}
	}
	return ends;
}

/*
Given an angularly ordered tree (vertices), construct a set of contours where
each contour describe a branch of the tree.
*/
vector<ContourEQW> 
extractPathTree(vector<Vertex<FragmentInfo>*>& vertices, 
				vector<ContourEQW>& contours)
{
	vector<ContourEQW> tree;
	vector<Node<Vertex<FragmentInfo>*>*> vnodes;
	for(int i=0; i<vertices.size(); ++i)
	{
		vnodes.push_back(makeset(vertices[i]));
	}

	vector<Vertex<FragmentInfo>*> ends = findEndPoints(vertices);
	for(int i=0; i<ends.size(); ++i)
	{
		Vertex<FragmentInfo>* p = ends[i];
		vector<Vertex<FragmentInfo>*> path;
		while(true)
		{
			path.push_back(p);
			if(p->pi == NULL)
			{
				break;
			}
			int k = distance(vertices.begin(), find(vertices.begin(), vertices.end(), p->pi));
			if(k < vnodes.size() && vnodes[k]->parent == vnodes[k])
			{
				merge(vnodes[i], vnodes[k]);
				p = p->pi;
			}
			else
			{
				path.push_back(p->pi);
				break;
			}
		}
		tree.push_back(stitchContour(path, contours));
	}

	for(int i=0; i<vnodes.size(); ++i)
	{
		delete vnodes[i];
		vnodes[i] = 0;
	}
	return tree;
}

