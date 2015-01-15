#ifndef _DISJOINT_SET_H_
#define _DISJOINT_SET_H_

template<class T>
struct Node
{
	Node(T k)
	{
		key=k;
		rank=0;
		parent = this;
	}
	void Reset()
	{
		rank=0; 
		parent = this;
	}
	T key;
	int rank;
	Node* parent;
};

template<class T>
Node<T>* makeset(T t)
{
	return new Node<T>(t);
}

template<class T>
Node<T>* merge(Node<T>* p, Node<T>* q)
{
	Node<T>* rp = findset(p);
	Node<T>* rq = findset(q);
	if(rp == rq)
	{
		return rp;
	}
	else if(rp->rank < rq->rank)
	{
		rp->parent = rq;
		return rp;
	}
	else if(rp->rank > rq->rank)
	{
		rq->parent = rp;
		return rq;
	}
	else 
	{
		rq->parent = rp;
		rp->rank++;
		return rq;
	}
}

template<class T>
Node<T>* 
findset(Node<T>* p)
{
	if(p != p->parent)
	{
		p->parent = findset(p->parent);
	}
	return p->parent;
}

#include <vector>
#include <algorithm>
using namespace std;

template<class T>
std::vector<Node<T>*>
clusters(const std::vector<Node<T>*>& vnodes)
{
	std::vector<Node<T>*> unique;
	for(int i=0; i<vnodes.size(); ++i)
	{
		Node<T>* r = findset(vnodes[i]);
		if(find(unique.begin(), unique.end(), r) == unique.end())
		{
			unique.push_back(r);
		}
	}

	return unique;
}

#include <cassert>
template<class T>
bool
removeNode(Node<T>* node, std::vector<Node<T>*>& vnodes)
{
	std::vector<Node<T>*>::iterator it = find(vnodes.begin(), vnodes.end(), node);
	assert(it != vnodes.end()); // return false; //node is not in the vector

	std::vector<Node<T>*> cluster;
	for(int i=0; i<vnodes.size(); ++i)
	{
		if(vnodes[i] == node) continue;

		Node<T>* r = findset(vnodes[i]);
		if(r == node)
		{
			cluster.push_back(vnodes[i]);
		}
	}
	if(node->parent == node)
	{
		if(cluster.empty() == false)
		{
			for(int i=0; i<cluster.size(); ++i)
			{
				cluster[i]->parent = cluster[0];
			}
		}
	}
	else
	{
		if(cluster.empty() == false)
		{
			for(int i=0; i<cluster.size(); ++i)
			{
				cluster[i]->parent = node->parent;
			}
		}
	}
	delete *it;
	vnodes.erase(it);
	return true;
}

#endif /* _DISJOINT_SET_H_ */