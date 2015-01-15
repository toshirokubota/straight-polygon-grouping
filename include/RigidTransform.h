#ifndef ___RIGID_TRANSFORM_H___
#define ___RIGID_TRANSFORM_H___
#include <myDataType.h>
#include <vector>
struct RigidTransform
{
public:
	RigidTransform()
	{
		translation[0] = translation[1] = 0;
		rotation = 0;
		scale[0] = scale[1] = 0;
	}
	RigidTransform(real a, real b, real c, real d, real e, real f)
	{
		translation[0] = c; translation[1] = f;
	}
	RigidTransform(std::vector<real>& p): RigidTransform(p[0], p[1], p[2], p[3], p[4], p[5])
	{
	}

	real translation[2];
	real rotation;
	real scale[2];
};


#endif /* ___RIGID_TRANSFORM_H___ */