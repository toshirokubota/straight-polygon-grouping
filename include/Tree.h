#ifndef ___TREE_H___
#define ___TREE_H___

#include <Graph.h>

template <class T>
class TreeNode: public Vertex<T>
{
public:
	TreeNode(T key): Vertex<T>(key)
	{
		parent = NULL;
	}
	vector<TreeNode<T>*> children;
	TreeNode<T>* parent;
};

#endif /* ___TREE_H___ */

