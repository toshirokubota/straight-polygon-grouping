#ifndef ___PRIM_H___
#define ___PRIM_H___

#include <vector>
using namespace std;
#include <Graph.h>

template <class T>
vector<Edge<T>*>
Prim(Vertex<T>* u, vector<Vertex<T>*> vertices);

#include "..\src\Prim.cpp"

#endif /* ___PRIM_H___ */