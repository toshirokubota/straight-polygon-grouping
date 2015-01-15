// Prim.cpp : Defines the entry point for the console application.
//

#include <Prim.h>

template<class T>
Vertex<T>*
ExtractMinimum(vector<Vertex<T>*>& Q)
{
	if(Q.empty()) return NULL;
	vector<Vertex<T>*>::iterator it = Q.begin();
	for(vector<Vertex<T>*>::iterator i=Q.begin()+1; i!=Q.end(); ++i)
	{
		if((*i)->d < (*it)->d)
		{
			it = i;
		}
	}
	Vertex<T>* u = *it;
	Q.erase(it);
	return u;
}

template <class T>
vector<Edge<T>*>
Prim(Vertex<T>* u, vector<Vertex<T>*> vertices)
{
	vector<Edge<T>*> mst;
	vector<Vertex<T>*> Q;
	//Step 1: Initialize all vertices. Call Reset(). Then, add it to Q.
	for(int i=0; i<vertices.size(); ++i)
	{
		vertices[i]->Reset();
		Q.push_back(vertices[i]);
	}
	//Step 2: Set the d value of the source (u) to 0. 
	u->d = 0;

	while(Q.empty() == false) 
	{
		//Extract a vertex with the minimum d
		Vertex<T>* u = ExtractMinimum(Q);
		if(u->pi)
		{
			Edge<T>* e = u->findEdge(u->pi);
			assert(e);
			mst.push_back(e);
		}
		//Update each adjacent edge.
		for(int i=0; i<u->aList.size(); ++i)
		{
			Edge<T>* e = u->aList[i];
			Vertex<T>* v = e->v;
			if(v->color == White)
			{
				if(v->d > e->w)
				{
					v->d = e->w;
					v->pi = u;
				}
			}
		}
		u->color = Black;
	}
	return mst;
}


