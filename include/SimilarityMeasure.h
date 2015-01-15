#ifndef ___SIMILARITY_MEASURE_H___
#define ___SIMILARITY_MEASURE_H___

#include <ContourEQW.h>

ContourEQW Interpolate(const ContourEQW& a);

float
ComputeSimilarity(const ContourEQW& a, const ContourEQW& b);

#endif /* ___SIMILARITY_MEASURE_H___ */
