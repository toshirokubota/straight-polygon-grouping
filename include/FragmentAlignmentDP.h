#ifndef ___FRAGMENT_ALIGNMENT_DP_H___
#define ___FRAGMENT_ALIGNMENT_DP_H___
#include <ContourEQW.h>

/*
This computes an alignment cost with the orientation of the fragment into consideration.
*/
float
StructuralDesimilarityDP(const ContourEQW& tc, const ContourEQW& sc, float skip, float repeat);

/*
This computes an alignment cost without the orientation of the fragment into consideration.
Thus, it is a weaker form than what is computed by StructuralDesimilarityDP.
*/
float
RandomDesimilarityDP(const ContourEQW& tc, const ContourEQW& sc, float skip, float repeat);

#endif /* ___FRAGMENT_ALIGNMENT_DP_H___ */