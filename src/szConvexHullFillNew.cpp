#ifdef MEX_DLL
#include <mex.h>
#endif
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;
#include <stdio.h>
#include <szDefaultParam.h>
#include <szMiscOperations.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szGeometry.h>

/*
To speed up, only provide volume surface to the convex hull process.
*/
void
preConvexHull(vector<unsigned char>& S,   //OUTPUT - must be initialized to zero
              const vector<unsigned char>& C,   //INPUT - segmentation map
              const int* dims) {  
  int i,j,k;
  
  for(k=0; k<dims[2]; ++k) 
  {
    for(j=0; j<dims[1]; ++j) 
    {
      for(i=0; i<dims[0]; ++i) 
      {
        if(onSurface3(C,i,j,k,dims)) {
          {
            SetData3(S,i,j,k,dims[0],dims[1],dims[2],(unsigned char)1);
          }
        }
      }
    }
  }
}

/*
Segmentation is the intersection of convex-hull and the original foreground segmentation.
This procedure takes care of the taks.
In a rare case, the cluster is broken up into multiple pieces after the above operation.
The function then choose the largest connected component.
*/
void
postConvexHull(vector<unsigned char>& B,    //OUTPUT - must be initialized to zero
               const vector<unsigned char>& C,    //INPUT - convex-hull map
               const vector<unsigned char>& L,   //INPUT - foreground segmentation map
               const int* dims)             //volume size
{
  int i;
  for(i=0; i<B.size(); ++i) 
  {
    if(C[i])
	{
		if(L[i])
			B[i] = ForegroundColor;
		else
			B[i] = 0;
	}
    else
      B[i] = 0;
  }
}

/*
Make sure the current configuration is correct.
This is for debugging purpose.
*/
bool
IntegrityCheck(const vector<szFace*>& v_pFaces) {
  for(int i=0; i<v_pFaces.size(); ++i) {
    szEdge* peBegin = v_pFaces[i]->getOuterComponent();
    szEdge* pe=peBegin;
    do {
      try {
        szVertex* pv = pe->getOrigin();
        pe=pe->getNext();
      }
      catch(...) {
        return false;
      }
    }
    while(pe!=peBegin);
    do {
      try {
        szVertex* pv = pe->getOrigin();
        pe=pe->getPrev();
      }
      catch(...) {
        return false;
      }
    }
    while(pe!=peBegin);
    do {
      szFace* pf = pe->getIncidentFace();
      szVertex* pv2=pe->getTwin()->getNext()->getNext()->getOrigin();
      if(CoPlaner(pv2,pf)) {
        cout << "Two adjacent faces cannot be co-linear.\n";
        cout << "Something is terribly wrong.\n";
        return false;
      }
      szVertex* pv = pe->getOrigin();
      pe=pe->getPrev();
    }
    while(pe!=peBegin);
    szFace* pf=v_pFaces[i];
    peBegin = pf->getOuterComponent()->getTwin()->getIncidentFace()->getOuterComponent();
    pe=peBegin;
    do {
      try {
        szVertex* pv = pe->getOrigin();
        pe=pe->getNext();
      }
      catch(...) {
        return false;
      }
    }
    while(pe!=peBegin);
    do {
      try {
        szVertex* pv = pe->getOrigin();
        pe=pe->getPrev();
      }
      catch(...) {
        return false;
      }
    }
    while(pe!=peBegin);
  }
  return true;
}

/*
Given a convex hull, it performs a half-plane segmentation.
*/
void
HalfPlanePaint(vector<unsigned char>& S, //segmentation result (OUTPUT)
               const vector<szFace*>& v_pFaces,//faces that define the convex hull (INPUT)
               int ndim,                  //the number of dimensions. Must be 3 (INPUT)
               const int* dims            //the size across each dimension
               )	 {
  
  vector<int> vlowerb(3,MAX_INTEGER);
  vector<int> vupperb(3,MIN_INTEGER);
  int i,j,k,n;
  //first compute a bounding box of the convex hull
  for(i=0; i<v_pFaces.size(); ++i) {
    szEdge* pe=v_pFaces[i]->getOuterComponent();
    szEdge* peBegin=pe;
    do {
      szVertex* pv = pe->getOrigin();
      double x=pv->getX();
      double y=pv->getY();
      double z=pv->getZ();
      vlowerb[0] = Min((int)x,vlowerb[0]);
      vupperb[0] = Max((int)x,vupperb[0]);
      vlowerb[1] = Min((int)y,vlowerb[1]);
      vupperb[1] = Max((int)y,vupperb[1]);
      vlowerb[2] = Min((int)z,vlowerb[2]);
      vupperb[2] = Max((int)z,vupperb[2]);
      pe=pe->getNext();
    }
    while(pe!=peBegin);
  }

  double slack=0.01;
  memset(&(S[0]),0,S.size()); //initialize
  for(n=0; n<v_pFaces.size(); n++) {
    vector<int> vnorm = v_pFaces[n]->getNormal();
    szVertex* pv = v_pFaces[n]->getOuterComponent()->getOrigin();
    int x = pv->getX();
    int y = pv->getY();
    int z = pv->getZ();
    //printf("%d: Normal: %d %d %d\tOrigin: %d %d %d\n", n, vnorm[0],vnorm[1],vnorm[2],x,y,z);

    //only look within a bounding box
    for(i=vlowerb[2]; i<=vupperb[2]; ++i) {
      for(j=vlowerb[1]; j<=vupperb[1]; ++j) {
        for(k=vlowerb[0]; k<=vupperb[0]; ++k) {
          if(n==0) //initialize it to 1
            SetData3(S,k,j,i,dims[0],dims[1],dims[2],(unsigned char)1);
          else if(GetData3(S,k,j,i,dims[0],dims[1],dims[2],(unsigned char)0)==0)
            continue;

          //if the dot product is positive, then it is outside the convex hull
          double dp = vnorm[0]*(k-x)+vnorm[1]*(j-y)+vnorm[2]*(i-z);
          if(dp>slack)
            SetData3(S,k,j,i,dims[0],dims[1],dims[2],(unsigned char)0); //set it to zero once for all
        }
      }
    }
  }

  int count=0;
  for(i=0; i<S.size(); ++i) {
    if(S[i])
      count++;
  }
  //printf("There are %d non-zero elements after half-plane paint.\n", count);
}

void
RandomizeVertices(vector<szVertex*>& v) {
  int i;
  int nitems=v.size();
  for(i=0; i<v.size(); ++i) {
    int k=(rand() % nitems) + i;
    szVertex* ptemp = v[i];
    v[i]=v[k];
    v[k]=ptemp;
    nitems--;
  }
}

/*
This routine finds 4 vertices that are not co-planer.  The vertices will form an initial convex hull.
*/
bool
FindInitialHull(vector<szVertex*>& v_pVertices,
                vector<szFace*>& v_pFaces) {
  szVertex* p1=v_pVertices[0];
  szVertex* p2=v_pVertices[1];
  szVertex* p3=NULL;
  int i, j;
  for(i=2; i<v_pVertices.size(); ++i) {
    if(!CoLinear(p1,p2,v_pVertices[i])) {
      p3=v_pVertices[i];
      v_pVertices[i]=v_pVertices[2];
      v_pVertices[2]=p3; //move the 3rd point to the 3rd element of the array
      break;
    }
  }
  if(p3==NULL) {
    //printf("No non-coliner points found.\n");
    return false;
  }

  szVertex* p4=NULL;
  for(i=3; i<v_pVertices.size(); ++i) {
    if(!CoPlaner(p1,p2,p3,v_pVertices[i])) {
      p4=v_pVertices[i];
      v_pVertices[i]=v_pVertices[3];
      v_pVertices[3]=p4; //move the 4th point to the 4th element of the array
      break;
    }
  }
  if(p4==NULL) {
    //printf("No non-coplaner points found.\n");
    return false;
  }

  szFace* pf[4];
  if(getVolume(p1,p2,p3,p4)<0)
    pf[0] = MakeNewFace(p1,p2,p3);
  else
    pf[0] = MakeNewFace(p3,p2,p1);
  if(getVolume(p2,p3,p4,p1)<0)
    pf[1] = MakeNewFace(p2,p3,p4);
  else
    pf[1] = MakeNewFace(p4,p3,p2);
  if(getVolume(p3,p4,p1,p2)<0)
    pf[2] = MakeNewFace(p3,p4,p1);
  else
    pf[2] = MakeNewFace(p1,p4,p3);
  if(getVolume(p4,p1,p2,p3)<0)
    pf[3] = MakeNewFace(p4,p1,p2);
  else
    pf[3] = MakeNewFace(p2,p1,p4);
  v_pFaces.push_back(pf[0]);
  v_pFaces.push_back(pf[1]);
  v_pFaces.push_back(pf[2]);
  v_pFaces.push_back(pf[3]);

  for(i=0; i<v_pFaces.size(); ++i) {
    szEdge* pe = v_pFaces[i]->getOuterComponent();
    if(pe==NULL) {
      //we should not enter this as all faces should have valid outer components
      pe = v_pFaces[i]->getInnerComponent();
      if(pe==NULL)
        continue; //this should never happen...
    }
    szEdge* peBegin = pe;
    do {
      if(pe->getTwin()==NULL) {
        for(j=(i+1)%v_pFaces.size(); j!=i; j=(j+1)%v_pFaces.size()) {
          szEdge* peTwin = FindTwin(v_pFaces[j],pe);
          if(peTwin!=NULL) {
            pe->setTwin(peTwin);
            break;
          }
        }
        if(j==i) {
          cout << "Could not find twin for an edge " << pe << endl;
          cout << "Something is terribly wrong!" << endl;
        }
      }
      pe = pe->getNext();
    }
    while(pe!=peBegin);
  }
  
  return true;
}

/*
Convlict graph contains a pair of vertex and face.  A face is conflicting to a vertex if addition of the vertex
to the convex hull will remove the face from the convex hull.
Initially, there are four faces.  This routine walk through each vertex that are not a part of the initial four
vertices, and add those that are outside the convex hull and associate each vertex with faces that are 
visible from the vertex.
*/
void
InitializeConflictGraph(szConflictGraph& graph, //conflict graph
                        const vector<szVertex*>& v_pVertices, //the set of vertices - first 4 formes the initla hull
                        const vector<szFace*>& v_pFaces //a set of four faces in the initial hull
                        ) {

  graph.m_lpFace.clear();
  graph.m_lpVertex.clear();

  int i,j;
  for(i=4; i<v_pVertices.size(); ++i) {
    for(j=0; j<v_pFaces.size(); ++j) {
      if(IsVisible(v_pVertices[i],v_pFaces[j])) {
        graph.m_lpFace.push_back(v_pFaces[j]);
        graph.m_lpVertex.push_back(v_pVertices[i]);
        //cout << *v_pVertices[i] << " -> Face #" << j << endl;
      }
    }
  }
}

vector<szFace*>
FindConflictingFaces(const szConflictGraph& graph,
                     const szVertex* pNewVertex) {
  vector<szFace*> v_pConflictFaces;
  list<szVertex*>::const_iterator lv = graph.m_lpVertex.begin();
  list<szFace*>::const_iterator lf = graph.m_lpFace.begin();
  for(; lv!=graph.m_lpVertex.end(); ++lv,++lf) {
    if(*lv == pNewVertex) { //conflict found
      v_pConflictFaces.push_back(*lf);
    }
  }

  //cout << "There are " << v_pConflictFaces.size() << " conflict faces.\n";
  //for(int i=0; i<v_pConflictFaces.size(); ++i)
  //  cout << *(v_pConflictFaces[i]) << endl;
  return v_pConflictFaces;
}

/*
This routine will look through conflicting faces and form a set of horizontal edges.
Horizontal edge divide visible face and invisible face.
*/
vector<szEdge*>
FindHorizontalEdges(const szVertex* pNewVertex, 
                    const vector<szFace*>& v_pConflictFaces) {

  vector<szEdge*> v_pHorizontalEdges;
  vector<szFace*>::const_iterator vf = v_pConflictFaces.begin();
  for(; vf!=v_pConflictFaces.end(); vf++) {
    szEdge* pe=(*vf)->getOuterComponent();
    szEdge* firstEdge = pe;
    do {
      //check if the face that shares this edge is visible
      szFace* pf = pe->getTwin()->getIncidentFace();
      if(find(v_pConflictFaces.begin(),v_pConflictFaces.end(),pf)<v_pConflictFaces.end()) {
        //this face is also in the conflicting face list - thus this edge is not a horizontal one
      }
      else {
        v_pHorizontalEdges.push_back(pe);
      }
      pe = pe->getNext();
    }
    while (pe!=firstEdge);
  }

  /*szEdge* peBegin = v_pHorizontalEdges[0];
  szEdge* pe=peBegin;
  do {
    if(CoLinear(pe->getOrigin(),pe->getNext()->getOrigin(),pe->getNext()->getNext()->getOrigin())) {
      cout << "Horizontal edges cannot be co-linera." << endl;
      cout << "Something is terribly wrong!" << endl;
    }
    pe=pe->getNext();
  }
  while(pe!=peBegin);*/

  //cout << "There are " << v_pHorizontalEdges.size() << " horizontal edges.\n";
  //for(int i=0; i<v_pHorizontalEdges.size(); ++i)
  //  cout << v_pHorizontalEdges[i] << "\t";
  //cout << endl;

  return v_pHorizontalEdges;
}

/*
This routine removes any bi-partite entries in the conflict graph that are associated with the new vertex.
We then add new entries to the graph that are associated with the newly added faces.
Note: only vertices already present in the conflict graph can form new conflicts to the new faces.
*/
void
UpdateConflictGraph(szConflictGraph& graph, //conflict graph
                    const szVertex* pNewVertex,   //newly added vertex
                    const vector<szFace*>& v_pNewFaces, //newly added faces
                    const vector<szFace*>& v_pConflictFaces //faces that are conflicting to the new vertex
                    ) {
  //first remove those entries that correspond to faces being removed
  //I am having problem with erase, so I clear them first and rebuilding it
  list<szVertex*> lvertices = graph.m_lpVertex;
  list<szFace*> lfaces = graph.m_lpFace;
  graph.m_lpFace.clear();
  graph.m_lpVertex.clear();
  list<szVertex*>::iterator lv;
  list<szFace*>::iterator lf;
  for(lf=lfaces.begin(),lv=lvertices.begin(); lf!=lfaces.end(); ++lf, ++lv) { 
    if(find(v_pConflictFaces.begin(),v_pConflictFaces.end(),*lf)==v_pConflictFaces.end()) {
      graph.m_lpVertex.push_back(*lv);
      graph.m_lpFace.push_back(*lf);
    }
  }
  
  //next add entries associated with newly created faces
  //first collect those vertices that are in the conflict graph
  vector<szVertex*> v_pUniqueVertices;
  for(lv = lvertices.begin(); lv!=lvertices.end(); ++lv) {
    if(find(v_pUniqueVertices.begin(),v_pUniqueVertices.end(),*lv)==v_pUniqueVertices.end())
      v_pUniqueVertices.push_back(*lv);
  }

  vector<szFace*>::const_iterator vf;
  for(vf = v_pNewFaces.begin(); vf<v_pNewFaces.end(); ++vf) {
    //we only need to look at those vertices already in the conflict graph
    vector<szVertex*>::const_iterator vv;
    for(vv = v_pUniqueVertices.begin(); vv!=v_pUniqueVertices.end(); ++vv) { 
      if(IsVisible(*vv,*vf)) {
        //cout << "Face [" << **vf << "] is visible from " << **vv << endl;
        graph.m_lpVertex.push_back(*vv);
        graph.m_lpFace.push_back(*vf);
      }
      //else 
        //cout << "Face [" << **vf << "] is NOT visible from " << **vv << endl;
    }
  }
}

/*
Given a new vertex and a set of horizontal edges, this routine will construct new faces.
There are some occasions where two newly formed faces are co-planer and adjacent to each other.
In such case, we merge the two faces into one.
*/
vector<szFace*>
MakeNewVisibleFaces(szVertex* pv, //new vertex
                    const vector<szEdge*>& v_pHorizontalEdges //a set of horizontal edges
                    ) {
  int i;
  int numE = v_pHorizontalEdges.size();
  //reorder horizontal edges so that they form a cycle
  szVertex* pvBegin=v_pHorizontalEdges[0]->getOrigin();
  szVertex* pvt = v_pHorizontalEdges[0]->getNext()->getOrigin();
  vector<szEdge*> v_pReordered;
  v_pReordered.push_back(v_pHorizontalEdges[0]);
  do {
    vector<szEdge*>::const_iterator pp;
    for(pp=v_pHorizontalEdges.begin(); pp<v_pHorizontalEdges.end(); ++pp) {
      if(pvt == (*pp)->getOrigin()) {
        v_pReordered.push_back(*pp);
        pvt=(*pp)->getNext()->getOrigin();
        break;
      }
    }
    if(pp == v_pHorizontalEdges.end()) {
      cout << "Cannot form a cycle out of horizontal edges..." << endl;
      cout << "There is something terribly wrong!\n";
      vector<szFace*> v;
      return v; //return an empty vector.
      //this should not happen;
    }
  }
  while(pvt!=pvBegin);

  /*szEdge* tmparray[10];
  for(i=0; i<v_pReordered.size() && i<10; ++i) {
    tmparray[i]=v_pReordered[i];
  }*/

  vector<szFace*> v_pFaces;
  for(i=0; i<v_pReordered.size(); ++i) {
    szFace* pf = MakeNewFace(pv,v_pReordered[i]);
    if(i>0 && CoPlaner(v_pReordered[i]->getPrev()->getOrigin(),pf)) {
      MergeFaces(*(v_pFaces.end()-1),pf);
      delete pf;
    }
    else
      v_pFaces.push_back(pf);
  }
  if(CoPlaner((*(v_pReordered.end()-1))->getOrigin(),*(v_pFaces.begin()))) {
    MergeFaces(*(v_pFaces.begin()),*(v_pFaces.end()-1));
    delete v_pFaces[v_pFaces.size()-1];
    v_pFaces.erase(v_pFaces.end()-1);
  }

  return v_pFaces;
}

/*
Now we know which faces to be removed and new faces to be added.
We need to do the following to maintain the integrity of the convex hull.
1. check if some new faces are co-planer to existing ones (that are adjacent to horizontal edges).
   if so, then merge them
2. check if two edges in a newly merged face are co-linear.  If so, merge them.
3. assign twin edges to new edges.
4. reassign twin edges for edges that were adjacent to horizontal edges (the horizontal edges will be removed)
5. remove those conflicting faces
*/
bool
UpdateGeometricalPrimitives(szVertex* pNewVertex,
                            vector<szFace*>& v_pNewFaces,
                            const vector<szFace*>& v_pFace2Remove,
                            const vector<szVertex*>& v_pVertices,
                            vector<szFace*>& v_pFaces,
                            bool trackOn) {
  int i,j;

  //check if some new faces are co-planer to existing ones
  //by construction of MakeNewFace(szVertex*, szEdge*), the OuterComponent of a new face is adjacent to
  //an existing face. This is the only one that can be co-planer to the new face
  vector<szFace*> v_pUpdatedFaces;
  //ofstream log("C:\\log.txt",ios::app);
  //log << *pNewVertex << endl;
  for(i=0; i<v_pNewFaces.size(); ++i) {
    //log << *(v_pNewFaces[i]) << endl;
    szFace* pf1=v_pNewFaces[i];
    szEdge* pe1=pf1->getOuterComponent();
    szEdge* pe2=pe1->getTwin();//twin is already set for this edge in MakeNewFaces
    szFace* pf2=pe2->getIncidentFace();
    if(CoPlaner(pNewVertex,pf2)) { 
      //log << "The new surface is co-planer to an existing one.\n";
      //merge to the existing one
      MergeFaces(pf2,pf1);
      //we no longer need the new face - replace with the old one that is just being updated
      delete pf1;
      v_pNewFaces[i]=NULL; //mark it as being deleted
      
      //add this modified surface into an array if it has not been done.
      if(find(v_pUpdatedFaces.begin(),v_pUpdatedFaces.end(),pf2)==v_pUpdatedFaces.end())
        v_pUpdatedFaces.push_back(pf2);
    }
    else {
      v_pFaces.push_back(pf1); //keep the faces that are on the convex hull
      v_pUpdatedFaces.push_back(pf1);
      pe2->setTwin(pe1); //reassign a twin edge for the edge adjacent to the new surface
      //log << "The new surface is NOT co-planer to an existing one.\n";
    }
  }

  //clean up deleted new faces and associated horizontal edges
  for(i=v_pNewFaces.size()-1; i>=0; --i) {
    if(v_pNewFaces[i]==NULL) {
      v_pNewFaces.erase(v_pNewFaces.begin()+i);
    }
  }

  //for all updated faces, remove colinear edges
  for(i=0; i<v_pUpdatedFaces.size(); ++i) {
    RemoveColinearEdges(v_pUpdatedFaces[i]);
  }

  //Now assign a twin edge to those who need it
  vector<szEdge*> v_pNoTwins;
  //to speed up search, put all edges who need twin into an array
  int count=0;
  //cout << "Checking who lacks twin in " << v_pUpdatedFaces.size() << " updated surfaces.\n";
  for(i=0; i<v_pUpdatedFaces.size(); ++i) {
    szFace* pf=v_pUpdatedFaces[i];
    szEdge* pe=pf->getOuterComponent();
    szEdge* peBegin=pe;
    do {
      if(pe->getTwin()==NULL) {
        v_pNoTwins.push_back(pe);
        //log << "NoTwin: " << pe << ": " << *pe << endl;
        count++;
      }

      pe=pe->getNext();
    }
    while(pe!=peBegin);
  }
  //cout << "There are " << count << "(" << v_pNoTwins.size() << ") edges that need twin.\n";
  for(i=0; i<count-1; ++i) {
    szEdge* pe=v_pNoTwins[i];
    if(pe->getTwin()==NULL) {
      for(j=i+1; j<count; j++) {
        szEdge* pe2=v_pNoTwins[j];
        if(IsTwin(pe,pe2)) {
          pe->setTwin(pe2);
          pe2->setTwin(pe);
          //cout << "Twin assignment: " << i << " - " << j << endl;
          break;
        }
      }
      if(j==v_pNoTwins.size()) {
        //log << "Could not find twin for an edge " << pe << endl;
        //log << "Something is terribly wrong!" << endl;
        /*for(int ii=0; ii<v_pFaces.size(); ++ii) {
          szFace* pf2=v_pFaces[ii];
          szEdge* pe2=pf2->getOuterComponent();
          szEdge* peBegin=pe2;
          do {
            log << pe2 << ":" << *pe2 << endl;
            if(IsTwin(pe,pe2)) {
              log << "Twin found in!" << endl;
            }
            pe2=pe2->getNext();
          }
          while(pe2!=peBegin);
        }*/
        return false;
      }
    }
  }

  //finally remove those conflicting faces
  for(i=0; i<v_pFace2Remove.size(); ++i) {
    vector<szFace*>::iterator pf = find(v_pFaces.begin(),v_pFaces.end(),v_pFace2Remove[i]);
    if(pf == v_pFaces.end()) {
      //log << "Could not find a conflicting face in v_pFaces" << endl;
      //log << "Something is terribly wrong!" << endl;
      return false;
    }
    else 
      v_pFaces.erase(pf);
    delete v_pFace2Remove[i];
  }

  //log.close();
  return true;
}

/*
This is the meat of the ConvexHull routine.  The result is retured as a half-plane segmentation volume.
The implementation is applicable only to 3D volume data and is based on the book:
M. de Berg, M. van Kreveld, M. Overmars, and O. Schwarzkopf,
"Computational Geometry - Algorithms and Applications", Springer.
The algorithm is described in page 236-239.

OUTPUT:
S: result as a half-plane segmentation
INPUT:
L: non-convex segmentation volume
ndim: # of dimensions - must be 3
dims: dimensions in each of the 3 dimensions.
*/

bool
ConvexHull3D (vector<unsigned char>& S, 
              const vector<unsigned char>& L, 
              int ndim, 
              const int* dims)	 {

  vector<szVertex*> v_pVertices; 
  vector<szFace*> v_pFaces;
  szConflictGraph CG;

  int nump = 0;
  int i;
  for(i=0; i<L.size(); ++i) {
    if(L[i]) {
      vector<int> vsub = Ind2Sub(i,ndim,dims);
      szVertex* pv = new szVertex;
      pv->setX(vsub[0]);
      pv->setY(vsub[1]);
      pv->setZ(vsub[2]);
      v_pVertices.push_back(pv);
      nump++;
    }
  }
  if(nump < 4)
  {
	  printf("ConvexHull3D: At least 4 positions are required.  Only %d are found.\n", nump);
	  return false;
  }
  //printf("ConvexHull3D: There are %d positions in the volume.\n", nump);

  RandomizeVertices(v_pVertices);

  bool bFlag = FindInitialHull(v_pVertices,v_pFaces);
  if(bFlag == false) {
    printf("No non-coplaner vertices are found in %d voxels ...  No further processing is necessary for ConvexHull3D.\n", 
      v_pVertices.size());
    return false;
  }
  int numVerticesInHull = 4;
  InitializeConflictGraph(CG,v_pVertices,v_pFaces);

  /*cout << "Initial Convex Hull:\n";
  for(i=0; i<v_pFaces.size(); ++i) {
    cout << *v_pFaces[i] << endl;
  }*/

  bool bSuccess = true;
  bool track=false;
  for(i=numVerticesInHull; i<nump && bSuccess; ++i) {
    szVertex* pv=*(v_pVertices.begin()+i);
    //cout << "\n\n" << i << ": Vertex: " << *pv << endl;
    vector<szFace*> v_pConflictFaces=FindConflictingFaces(CG,pv);
    if(!v_pConflictFaces.empty()) {//there are conflicts
      vector<szEdge*> v_pHorizontalEdges=FindHorizontalEdges(pv,v_pConflictFaces);
      vector<szFace*> v_pNewFaces = MakeNewVisibleFaces(pv,v_pHorizontalEdges);
      bSuccess=UpdateGeometricalPrimitives(pv,v_pNewFaces,v_pConflictFaces,v_pVertices,v_pFaces,track);
      UpdateConflictGraph(CG,pv,v_pNewFaces,v_pConflictFaces);
    }
    else {
      //the point is inside the convex hull
    }
    //if(bSuccess)
    //  bSuccess=IntegrityCheck(v_pFaces);
  }

  /*cout << "Final Convex Hull:\n";
  for(i=0; i<v_pFaces.size(); ++i) {
    cout << *v_pFaces[i] << endl;
  }*/

  static int failedVolumes=1;
  if(!bSuccess || v_pFaces.size()<4) {
    printf("Convex Hull failed. %d %d\n", bSuccess, v_pFaces.size());
    //char filename[128];
    //sprintf(filename,"C:\\dumped%3.3d.data", failedVolumes);
    //FILE* fp=fopen(filename,"wb");
    //fwrite(L.begin(),L.size(),sizeof(unsigned char),fp);
    failedVolumes++;
  }
  else {
    //printf("Performing half-plane segmentation.\n");
    HalfPlanePaint(S,v_pFaces,ndim,dims);
  }

  for(i=0; i<v_pFaces.size(); ++i)
    delete v_pFaces[i];
  for(i=0; i<v_pVertices.size(); ++i)
    delete v_pVertices[i];

  return true;
}

