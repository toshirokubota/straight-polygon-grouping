//#include <szConnectedComponentC.h>
//#include <szMyNeighborOp.h>
#include <algorithm>
using namespace std;

template<class Item>
int
ConnectedComponentAnalysisBigger(vector<int>& C,  //has to be initialized to zero
                                 const vector<Item>& V,
                                 const vector<int>& nbh,
								 vector<int>& vcount,
                                 Item thres, 
                                 int ndim, 
                                 const int* dims) {
  int i,j;
  vector<int> lut(1,0); //look-up table
  int next_c=1;
  int nvoxels = numberOfElements(ndim,dims);
  if(C.size()<nvoxels || V.size()<nvoxels) {
    //mexPrintf("An input image does not contain enough data.\n");
    return 0;
  }

  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoords;
  for(i=0; i<nbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
    vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
  }

  for(i=0; i<nvoxels; ++i) {
    Item v=V[i];
    if (v>thres) {
      vector<int> vlabels;
      vector<int> vsub = Ind2Sub(i,ndim,dims);
      int count=0;
      int minL=nvoxels; //loose upper bound of labels
      for(j=0; j<nbh.size(); ++j) {
        if(!NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims))
          continue;
        int k=i+nbh[j];
        int lp = lut[C[k]];
        if(lp) {
          minL=Min(minL,lp);
          count++;
        }
        if(find(vlabels.begin(),vlabels.end(),lp)>=vlabels.end())
          vlabels.push_back(lp);
      }
      if(count) {
        C[i]=minL;
        //check for multiple labels merging at this point
        //if so, relabel them with the minimum number
        for(j=0; j<vlabels.size(); ++j) {
          if(vlabels[j]>minL) {
            for(int l=1; l<lut.size(); ++l) 
              if(lut[l]==vlabels[j])
                lut[l]=minL;
          }
        }
      }
      else {
        C[i]=next_c;
        lut.push_back(next_c);
        next_c++;
      }
    }
  }

  if(next_c==1) //no cluster found
    return 0;

  int max_label=0;
  for(i=1; i<lut.size(); ++i) {
    max_label=Max(max_label,lut[i]);
  }

  //reassign labels in consecutive order
  next_c=1;
  for(i=1; i<=max_label; ++i) {
    vector<int>::iterator iter = lut.begin();
    bool found=false;
    do {
      iter=find(iter, lut.end(), i);
      if(iter<lut.end()) {
        *iter = next_c;
        iter++; //move one to search from the next item
        found = true;
      }
    }
    while(iter<lut.end());
    if(found)
      next_c++;
  }

  //relabel the data and count the number pixels in each component
  vcount = vector<int>(next_c+1, 0);
  for(i=0; i<nvoxels; ++i) {
    if(C[i]) {
      C[i] = lut[C[i]];
    }
	vcount[C[i]]++;
  }
  
  return next_c-1; //subtract 1 to disregard the background
}

template<class Item>
int
ConnectedComponentAnalysisSmaller(vector<int>& C,  //has to be initialized to zero
                                 const vector<Item>& V,
                                 const vector<int>& nbh,
								 vector<int>& vcount,
                                 Item thres, 
                                 int ndim, 
                                 const int* dims) {
  int i,j;
  vector<int> lut(1,0); //look-up table
  int next_c=1;
  int nvoxels = numberOfElements(ndim,dims);
  if(C.size()<nvoxels || V.size()<nvoxels) {
    //mexPrintf("An input image does not contain enough data.\n");
    return 0;
  }

  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoords;
  for(i=0; i<nbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
    vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
  }

  for(i=0; i<nvoxels; ++i) {
    Item v=V[i];
    if (v<thres) {
      vector<int> vlabels;
      vector<int> vsub = Ind2Sub(i,ndim,dims);
      int count=0;
      int minL=nvoxels; //loose upper bound of labels
      for(j=0; j<nbh.size(); ++j) {
        if(!NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims))
          continue;
        int k=i+nbh[j];
        int lp = lut[C[k]];
        if(lp) {
          minL=Min(minL,lp);
          count++;
        }
        if(find(vlabels.begin(),vlabels.end(),lp)>=vlabels.end())
          vlabels.push_back(lp);
      }
      if(count) {
        C[i]=minL;
        //check for multiple labels merging at this point
        //if so, relabel them with the minimum number
        for(j=0; j<vlabels.size(); ++j) {
          if(vlabels[j]>minL) {
            for(int l=1; l<lut.size(); ++l) 
              if(lut[l]==vlabels[j])
                lut[l]=minL;
          }
        }
      }
      else {
        C[i]=next_c;
        lut.push_back(next_c);
        next_c++;
      }
    }
  }

  if(next_c==1) //no cluster found
    return 0;

  int max_label=0;
  for(i=1; i<lut.size(); ++i) {
    max_label=Max(max_label,lut[i]);
  }

  //reassign labels in consecutive order
  next_c=1;
  for(i=1; i<=max_label; ++i) {
    vector<int>::iterator iter = lut.begin();
    bool found=false;
    do {
      iter=find(iter, lut.end(), i);
      if(iter<lut.end()) {
        *iter = next_c;
        iter++; //move one to search from the next item
        found = true;
      }
    }
    while(iter<lut.end());
    if(found)
      next_c++;
  }

  //relabel the data and count the number pixels in each component
  vcount = vector<int>(next_c+1, 0);
  for(i=0; i<nvoxels; ++i) {
    if(C[i]) {
      C[i] = lut[C[i]];
    }
	vcount[C[i]]++;
  }
  
  return next_c-1; //subtract 1 to disregard the background
}

template<class Item>
int
ConnectedComponentAnalysisEqual(vector<int>& C,  //has to be initialized to zero
                                const vector<Item>& V, 
                                const vector<int>& nbh,
								vector<int>& vcount,
                                Item value, 
                                int ndim, 
                                const int* dims) {
  int i,j;
  vector<int> lut(1,0); //look-up table
  int next_c=1;
  int nvoxels = numberOfElements(ndim,dims);
  if(C.size()<nvoxels || V.size()<nvoxels) {
    //mexPrintf("An input image does not contain enough data.\n");
    return 0;
  }

  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoords;
  for(i=0; i<nbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
    vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
  }

  for(i=0; i<nvoxels; ++i) {
    Item v=V[i];
    if (v==value) {
      vector<int> vlabels;
      vector<int> vsub = Ind2Sub(i,ndim,dims);
      int count=0;
      int minL=nvoxels; //loose upper bound of labels
      for(j=0; j<nbh.size(); ++j) {
        if(!NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims))
          continue;
        int k=i+nbh[j];
        int lp = lut[C[k]];
        if(lp) {
          minL=Min(minL,lp);
          count++;
        }
        if(find(vlabels.begin(),vlabels.end(),lp)>=vlabels.end())
          vlabels.push_back(lp);
      }
      if(count) {
        C[i]=minL;
        //check for multiple labels merging at this point
        //if so, relabel them with the minimum number
        for(j=0; j<vlabels.size(); ++j) {
          if(vlabels[j]>minL) {
            for(int l=1; l<lut.size(); ++l) 
              if(lut[l]==vlabels[j])
                lut[l]=minL;
          }
        }
      }
      else {
        C[i]=next_c;
        lut.push_back(next_c);
        next_c++;
      }
    }
  }

  if(next_c==1) //no cluster found
    return 0;

  int max_label=0;
  for(i=1; i<lut.size(); ++i) {
    max_label=Max(max_label,lut[i]);
  }

  //reassign labels in consecutive order
  next_c=1;
  for(i=1; i<=max_label; ++i) {
    vector<int>::iterator iter = lut.begin();
    bool found=false;
    do {
      iter=find(iter, lut.end(), i);
      if(iter<lut.end()) {
        *iter = next_c;
        iter++; //move one to search from the next item
        found = true;
      }
    }
    while(iter<lut.end());
    if(found)
      next_c++;
  }

  //relabel the data and count the number pixels in each component
  vcount = vector<int>(next_c+1, 0);
  for(i=0; i<nvoxels; ++i) {
    if(C[i]) {
      C[i] = lut[C[i]];
    }
	vcount[C[i]]++;
  }
  
  return next_c-1; //subtract 1 to disregard the background
}

template<class Item>
int
ConnectedComponentAnalysisEqualSizeThreshold(vector<int>& C,  //has to be initialized to zero
                                             const vector<Item>& V, 
                                             const vector<int>& nbh,
											 vector<int>& vcount,
                                             Item value, 
                                             int min_size,
                                             int max_size,
                                             int ndim, 
                                             const int* dims) {
  int i,j;
  vector<int> lut(1,0); //look-up table
  vector<int> vsize(1,0); //hold the cluster sizes 
  int next_c=1;
  int nvoxels = numberOfElements(ndim,dims);
  if(C.size()<nvoxels || V.size()<nvoxels) {
    //mexPrintf("An input image does not contain enough data.\n");
    return 0;
  }

  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoords;
  for(i=0; i<nbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
    vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
  }

  for(i=0; i<nvoxels; ++i) {
    Item v=V[i];
    if (v==value) {
      vector<int> vlabels;
      vector<int> vsub = Ind2Sub(i,ndim,dims);
      int count=0;
      int minL=nvoxels; //loose upper bound of labels
      for(j=0; j<nbh.size(); ++j) {
        if(!NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims))
          continue;
        int k=i+nbh[j];
        int lp = lut[C[k]];
        if(lp) {
          minL=Min(minL,lp);
          count++;
        }
        if(find(vlabels.begin(),vlabels.end(),lp)>=vlabels.end())
          vlabels.push_back(lp);
      }
      if(count) {
        C[i]=minL;
        //check for multiple labels merging at this point
        //if so, relabel them with the minimum number
        for(j=0; j<vlabels.size(); ++j) {
          if(vlabels[j]>minL) { //labels need to be merged
            for(int l=1; l<lut.size(); ++l) {
              if(lut[l]==vlabels[j]) {
                lut[l]=minL;
                if(vsize[l]) {
                  vsize[minL]+=vsize[l];
                  vsize[l]=0;
                }
              }
            }
          }
        }
      }
      else {
        C[i]=next_c;
        lut.push_back(next_c);
        vsize.push_back(1);
        next_c++;
      }
    }
  }

  if(next_c==1) //no cluster found
    return 0;

  //remove clusters whose size does not meet the requirement
  for(i=1; i<lut.size(); ++i) {
    if(vsize[i]<min_size || vsize[i]>max_size) {
      for(j=i; j<lut.size(); ++j) {
        if(lut[j]==i) {
          lut[j]=0; // set it as background
        }
      }
      vsize[i]=0;
    }
  }

  int max_label=0;
  for(i=1; i<lut.size(); ++i) {
    max_label=Max(max_label,lut[i]);
  }

  //reassign labels in consecutive order
  next_c=1;
  for(i=1; i<=max_label; ++i) {
    vector<int>::iterator iter = lut.begin();
    bool found=false;
    do {
      iter=find(iter, lut.end(), i);
      if(iter<lut.end()) {
        *iter = next_c;
        iter++; //move one to search from the next item
        found = true;
      }
    }
    while(iter<lut.end());
    if(found)
      next_c++;
  }

  //relabel the data and count the number pixels in each component
  vcount = vector<int>(next_c+1, 0);
  for(i=0; i<nvoxels; ++i) {
    if(C[i]) {
      C[i] = lut[C[i]];
    }
	vcount[C[i]]++;
  }  

  return next_c-1; //subtract 1 to disregard the background
}

template<class Item>
int
ConnectedComponentAnalysisEqual(vector<int>& C,  //has to be initialized to zero
                                const vector<Item>& V, 
								vector<int>& vcount,
                                int mode,
                                Item value, 
                                int ndim, 
                                const int* dims) {

  vector<int> nbh;
  int iRet;
  switch(mode) {
  case NeighborhoodFour: //four neighborhood generalized to N-dim
    nbh = MakeCausalFourNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisEqual(C,V,nbh,vcount,value,ndim,dims);
    break;
  case NeighborhoodEight: //eight neighborhood generalized to N-dim
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisEqual(C,V,nbh,vcount,value,ndim,dims);
    break;
  default:
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisEqual(C,V,nbh,vcount,value,ndim,dims);
    break;
  }

  return iRet;
}

template<class Item>
int
ConnectedComponentAnalysisEqualSizeThreshold(vector<int>& C,  //has to be initialized to zero
                                             const vector<Item>& V, 
											 vector<int>& vcount,
                                             int mode,
                                             Item value,
                                             int min_size,
                                             int max_size,
                                             int ndim, 
                                             const int* dims) {
  
  vector<int> nbh;
  int iRet;
  switch(mode) {
  case NeighborhoodFour: //four neighborhood generalized to N-dim
    nbh = MakeCausalFourNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisEqualSizeThreshold(C,V,nbh,vcount,value,min_size,max_size,ndim,dims);
    break;
  case NeighborhoodEight: //eight neighborhood generalized to N-dim
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisEqualSizeThreshold(C,V,nbh,vcount,value,min_size,max_size,ndim,dims);
    break;
  default:
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisEqualSizeThreshold(C,V,nbh,vcount,value,min_size,max_size,ndim,dims);
    break;
  }
  
  return iRet;
}

template<class Item>
int
ConnectedComponentAnalysisBigger(vector<int>& C,  //has to be initialized to zero
                                const vector<Item>& V, 
								vector<int>& vcount,
                                int mode,
                                Item value, 
                                int ndim, 
                                const int* dims) {

  vector<int> nbh;
  int iRet;
  switch(mode) {
  case NeighborhoodFour: //four neighborhood generalized to N-dim
    nbh = MakeCausalFourNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisBigger(C,V,nbh,vcount,value,ndim,dims);
    break;
  case NeighborhoodEight: //eight neighborhood generalized to N-dim
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisBigger(C,V,nbh,vcount,value,ndim,dims);
    break;
  default:
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisBigger(C,V,nbh,vcount,value,ndim,dims);
    break;
  }

  return iRet;
}

template<class Item>
int
ConnectedComponentAnalysisSmaller(vector<int>& C,  //has to be initialized to zero
                                const vector<Item>& V, 
								vector<int>& vcount,
                                int mode,
                                Item value, 
                                int ndim, 
                                const int* dims) {

  vector<int> nbh;
  int iRet;
  switch(mode) {
  case NeighborhoodFour: //four neighborhood generalized to N-dim
    nbh = MakeCausalFourNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisSmaller(C,V,nbh,vcount,value,ndim,dims);
    break;
  case NeighborhoodEight: //eight neighborhood generalized to N-dim
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisSmaller(C,V,nbh,vcount,value,ndim,dims);
    break;
  default:
    nbh = MakeCausalNeighborhood(ndim,dims);
    iRet = ConnectedComponentAnalysisSmaller(C,V,nbh,vcount,value,ndim,dims);
    break;
  }

  return iRet;
}
