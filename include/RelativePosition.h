#ifndef ___RELATIVE_POSITION_H___
#define ___RELATIVE_POSITION_H___
#include <ContourEQW.h>

enum RelativePosition {LeftSide, RightSide, FrontSide, BackSide, UndeterminedSide, UnknownSide};

RelativePosition
Side(float orientation, float direction);

RelativePosition
Side(ContourEQW& a, ContourEQW& b);

RelativePosition
FrontBack(ContourEQW& a, ContourEQW& b);

#endif /* ___RELATIVE_POSITION_H___ */