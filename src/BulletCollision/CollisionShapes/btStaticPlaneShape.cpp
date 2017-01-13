/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2009 Erwin Coumans  http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "btStaticPlaneShape.h"

#include "LinearMath/btTransformUtil.h"


btStaticPlaneShape::btStaticPlaneShape(const btVector3& planeNormal,btScalar planeConstant)
: btConcaveShape (), m_planeNormal(planeNormal.normalized()),
m_planeConstant(planeConstant),
m_localScaling(btScalar(1.),btScalar(1.),btScalar(1.))
{
	m_shapeType = STATIC_PLANE_PROXYTYPE;
	//	btAssert( btFuzzyZero(m_planeNormal.length() - btScalar(1.)) );
}


btStaticPlaneShape::~btStaticPlaneShape()
{
}


/// Zero out any entries that are approximately zero.
SIMD_FORCE_INLINE btVector3
btRoundZero(const btVector3 &vec, btScalar tol = 0.00001)
{
    btVector3 v = vec;
    for (int i = 0; i < 3; ++i)
    {
        if (btFabs(v[i]) <= tol)
            v[i] = 0;
    }

    return v;
}

/// Add a point to the bounding box.
SIMD_FORCE_INLINE void
btAddPoint(btVector3 &aabbMin, btVector3 &aabbMax, const btVector3 &p)
{
    aabbMin.setMin(p);
    aabbMax.setMax(p);
}
 
void btStaticPlaneShape::getAabb(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const
{
    const btVector3 center = t * (m_planeConstant * m_planeNormal);
    const btQuaternion r = shortestArcQuat(btVector3(0, 1, 0), t.getBasis() * m_planeNormal);
    const btVector3 v1 = btRoundZero(quatRotate(r, btVector3(1, 0, 0))) * BT_LARGE_FLOAT;
    const btVector3 v2 = btRoundZero(quatRotate(r, btVector3(0, 0, 1))) * BT_LARGE_FLOAT;
    const btVector3 v3 = btRoundZero(quatRotate(r, btVector3(0, 1, 0))) * BT_LARGE_FLOAT;

    aabbMin = aabbMax = center - v1;
    btAddPoint(aabbMin, aabbMax, center + v1);
    btAddPoint(aabbMin, aabbMax, center - v2);
    btAddPoint(aabbMin, aabbMax, center + v2);
    btAddPoint(aabbMin, aabbMax, center - v3);
}

void	btStaticPlaneShape::processAllTriangles(btTriangleCallback* callback,const btVector3& aabbMin,const btVector3& aabbMax) const
{

	btVector3 halfExtents = (aabbMax - aabbMin) * btScalar(0.5);
	btScalar radius = halfExtents.length();
	btVector3 center = (aabbMax + aabbMin) * btScalar(0.5);
	
	//this is where the triangles are generated, given AABB and plane equation (normal/constant)

	btVector3 tangentDir0,tangentDir1;

	//tangentDir0/tangentDir1 can be precalculated
	btPlaneSpace1(m_planeNormal,tangentDir0,tangentDir1);

	btVector3 supVertex0,supVertex1;

	btVector3 projectedCenter = center - (m_planeNormal.dot(center) - m_planeConstant)*m_planeNormal;
	
	btVector3 triangle[3];
	triangle[0] = projectedCenter + tangentDir0*radius + tangentDir1*radius;
	triangle[1] = projectedCenter + tangentDir0*radius - tangentDir1*radius;
	triangle[2] = projectedCenter - tangentDir0*radius - tangentDir1*radius;

	callback->processTriangle(triangle,0,0);

	triangle[0] = projectedCenter - tangentDir0*radius - tangentDir1*radius;
	triangle[1] = projectedCenter - tangentDir0*radius + tangentDir1*radius;
	triangle[2] = projectedCenter + tangentDir0*radius + tangentDir1*radius;

	callback->processTriangle(triangle,0,1);

}

void	btStaticPlaneShape::calculateLocalInertia(btScalar mass,btVector3& inertia) const
{
	(void)mass;

	//moving concave objects not supported
	
	inertia.setValue(btScalar(0.),btScalar(0.),btScalar(0.));
}

void	btStaticPlaneShape::setLocalScaling(const btVector3& scaling)
{
	m_localScaling = scaling;
}
const btVector3& btStaticPlaneShape::getLocalScaling() const
{
	return m_localScaling;
}
