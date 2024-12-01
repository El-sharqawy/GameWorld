#include "matrix.h"
#include <cmath>
#include <glm/gtx/quaternion.hpp>

void CMatrix4Df::InitRoatationX(const float fRotX, const bool bLeftRot)
{
	/* Rotate The Matrix in Left Handed Coordinate System Around X Axis */
	if (bLeftRot)
	{
		mat4[0][0] = 1.0f; mat4[0][1] = 0.0f; mat4[0][2] = 0.0f; mat4[0][3] = 0.0f;
		mat4[1][0] = 0.0f; mat4[1][1] = std::cosf(fRotX); mat4[1][2] = std::sinf(fRotX); mat4[1][3] = 0.0f;
		mat4[2][0] = 0.0f; mat4[2][1] = -std::sinf(fRotX); mat4[2][2] = std::cosf(fRotX); mat4[2][3] = 0.0f;
		mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
		return;
	}

	mat4[0][0] = 1.0f; mat4[0][1] = 0.0f; mat4[0][2] = 0.0f; mat4[0][3] = 0.0f;
	mat4[1][0] = 0.0f; mat4[1][1] = std::cosf(fRotX); mat4[1][2] = -std::sinf(fRotX); mat4[1][3] = 0.0f;
	mat4[2][0] = 0.0f; mat4[2][1] = std::sinf(fRotX); mat4[2][2] = std::cosf(fRotX); mat4[2][3] = 0.0f;
	mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
}

void CMatrix4Df::InitRoatationY(const float fRotY, const bool bLeftRot)
{
	/* Rotate The Matrix in Left Handed Coordinate System Around Y Axis */
	if (bLeftRot)
	{
		mat4[0][0] = std::cosf(fRotY); mat4[0][1] = 0.0f; mat4[0][2] = -std::sinf(fRotY); mat4[0][3] = 0.0f;
		mat4[1][0] = 0.0f; mat4[1][1] = 1.0f; mat4[1][2] = 0.0f; mat4[1][3] = 0.0f;
		mat4[2][0] = std::sinf(fRotY); mat4[2][1] = 0.0f; mat4[2][2] = std::cosf(fRotY); mat4[2][3] = 0.0f;
		mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
		return;
	}

	mat4[0][0] = std::cosf(fRotY); mat4[0][1] = 0.0f; mat4[0][2] = std::sinf(fRotY); mat4[0][3] = 0.0f;
	mat4[1][0] = 0.0f; mat4[1][1] = 1.0f; mat4[1][2] = 0.0f; mat4[1][3] = 0.0f;
	mat4[2][0] = -std::sinf(fRotY); mat4[2][1] = 0.0f; mat4[2][2] = std::cosf(fRotY); mat4[2][3] = 0.0f;
	mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
}

void CMatrix4Df::InitRoatationZ(const float fRotZ, const bool bLeftRot)
{
	/* Rotate The Matrix in Left Handed Coordinate System Around Z Axis */
	if (bLeftRot)
	{
		mat4[0][0] = std::cosf(fRotZ); mat4[0][1] = std::sinf(fRotZ); mat4[0][2] = 0.0f; mat4[0][3] = 0.0f;
		mat4[1][0] = -std::sinf(fRotZ); mat4[1][1] = std::cosf(fRotZ); mat4[1][2] = 0.0f; mat4[1][3] = 0.0f;
		mat4[2][0] = 0.0f; mat4[2][1] = 0.0f; mat4[2][2] = 1.0f; mat4[2][3] = 0.0f;
		mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
		return;
	}

	mat4[0][0] = std::cosf(fRotZ); mat4[0][1] = -std::sinf(fRotZ); mat4[0][2] = 0.0f; mat4[0][3] = 0.0f;
	mat4[1][0] = std::sinf(fRotZ); mat4[1][1] = std::cosf(fRotZ); mat4[1][2] = 0.0f; mat4[1][3] = 0.0f;
	mat4[2][0] = 0.0f; mat4[2][1] = 0.0f; mat4[2][2] = 1.0f; mat4[2][3] = 0.0f;
	mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
}

void CMatrix4Df::InitScaleTransform(const float fScaleX, const float fScaleY, const float fScaleZ)
{
	mat4[0][0] = fScaleX; mat4[0][1] = 0.0f; mat4[0][2] = 0.0f; mat4[0][3] = 0.0f;
	mat4[1][0] = 0.0f; mat4[1][1] = fScaleY; mat4[1][2] = 0.0f; mat4[1][3] = 0.0f;
	mat4[2][0] = 0.0f; mat4[2][1] = 0.0f; mat4[2][2] = fScaleZ; mat4[2][3] = 0.0f;
	mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
}

void CMatrix4Df::InitScaleTransform(const float fScale)
{
	InitScaleTransform(fScale, fScale, fScale);
}

void CMatrix4Df::InitScaleTransform(const SVector3Df& vScale)
{
	InitScaleTransform(vScale.x, vScale.y, vScale.z);
}

void CMatrix4Df::InitScaleTransform(const glm::vec3& vScale)
{
	InitScaleTransform(vScale.x, vScale.y, vScale.z);
}

/**
 * Initializes the matrix for a rotation transformation using the ZYX Euler rotation sequence.
 *
 * @param RotateX: The angle of rotation around the X-axis in degrees.
 * @param RotateY: The angle of rotation around the Y-axis in degrees.
 * @param RotateZ: The angle of rotation around the Z-axis in degrees.
 *
 * This function creates individual rotation matrices for the X, Y, and Z axes
 * and combines them in the order rz * ry * rx. This order applies the
 * rotations in the sequence: X-axis, Y-axis, then Z-axis, relative to the
 * local coordinate system.
 */
void CMatrix4Df::InitRotateTransform(const float fRotateX, const float fRotateY, const float fRotateZ)
{
	CMatrix4Df matX{}, matY{}, matZ{};

	const float x = ToRadian(fRotateX);
	const float y = ToRadian(fRotateY);
	const float z = ToRadian(fRotateZ);

	matX.InitRoatationX(x);
	matY.InitRoatationX(y);
	matZ.InitRoatationX(z);
	
	/* ZYX Euler rotation sequence */
	*this = matZ * matY * matX;
}

/**
 * Initializes the matrix for a rotation transformation using the XYZ Euler rotation sequence.
 *
 * @param RotateX: The angle of rotation around the X-axis in degrees.
 * @param RotateY: The angle of rotation around the Y-axis in degrees.
 * @param RotateZ: The angle of rotation around the Z-axis in degrees.
 *
 * This function creates individual rotation matrices for the X, Y, and Z axes
 * and combines them in the order rx * ry * rz. This order applies the
 * rotations in the sequence: Z-axis, Y-axis, then X-axis, relative to the
 * global coordinate system.
 */
void CMatrix4Df::InitRotateTransformZYX(const float fRotateX, const float fRotateY, const float fRotateZ)
{
	CMatrix4Df matX{}, matY{}, matZ{};

	const float x = ToRadian(fRotateX);
	const float y = ToRadian(fRotateY);
	const float z = ToRadian(fRotateZ);

	matX.InitRoatationX(x);
	matY.InitRoatationX(y);
	matZ.InitRoatationX(z);

	/* XYZ Euler rotation sequence */
	*this = matX * matY * matZ;
}

/**
 * Initializes the matrix for a rotation transformation using the ZYX Euler rotation sequence.
 *
 * @param vRotation: The 3D Vector holding rotation degrees.
 *
 * This function creates individual rotation matrices for the X, Y, and Z axes
 * and combines them in the order rz * ry * rx. This order applies the
 * rotations in the sequence: X-axis, Y-axis, then Z-axis, relative to the
 * local coordinate system.
 */
void CMatrix4Df::InitRotateTransform(const SVector3Df& vRotation)
{
	InitRotateTransform(vRotation.x, vRotation.y, vRotation.z);
}

/**
 * Initializes the matrix for a rotation transformation using the ZYX Euler rotation sequence.
 *
 * @param vRotation: The glm::vec3 holding rotation degrees.
 *
 * This function creates individual rotation matrices for the X, Y, and Z axes
 * and combines them in the order rz * ry * rx. This order applies the
 * rotations in the sequence: X-axis, Y-axis, then Z-axis, relative to the
 * local coordinate system.
 */
void CMatrix4Df::InitRotateTransform(const glm::vec3& vRotation)
{
	InitRotateTransform(vRotation.x, vRotation.y, vRotation.z);
}

void CMatrix4Df::InitRotateTransform(const SQuaternion& sQuat)
{
	const float yy2 = 2.0f * sQuat.y * sQuat.y;
	const float xy2 = 2.0f * sQuat.x * sQuat.y;
	const float xz2 = 2.0f * sQuat.x * sQuat.z;
	const float yz2 = 2.0f * sQuat.y * sQuat.z;
	const float zz2 = 2.0f * sQuat.z * sQuat.z;
	const float wz2 = 2.0f * sQuat.w * sQuat.z;
	const float wy2 = 2.0f * sQuat.w * sQuat.y;
	const float wx2 = 2.0f * sQuat.w * sQuat.x;
	const float xx2 = 2.0f * sQuat.x * sQuat.x;

	mat4[0][0] = -yy2 - zz2 + 1.0f;
	mat4[0][1] = xy2 + wz2;
	mat4[0][2] = xz2 - wy2;
	mat4[0][3] = 0.0f;

	mat4[1][0] = xy2 - wz2;
	mat4[1][1] = -xx2 - zz2 + 1.0f;
	mat4[1][2] = yz2 + wx2;

	mat4[2][0] = xz2 + wy2;
	mat4[2][1] = yz2 - wx2;
	mat4[2][2] = -xx2 - yy2 + 1.0f;

	mat4[3][0] = 0.0f;
	mat4[3][1] = 0.0f;
	mat4[3][2] = 0.0f;
	mat4[3][3] = 1.0f;
}

/**
 * Initializes the matrix for a rotation transformation using a quaternion.
 * @param sQuat: The quaternion representing the rotation to be applied from GLM.
 *
 * This function computes a 4x4 rotation matrix based on the given quaternion.
 * Quaternions provide a compact and efficient way to represent 3D rotations,
 * avoiding issues like gimbal lock that occur with Euler angles.
 *
 * The matrix is calculated directly from the quaternion's components:
 * - @quat.x: The X component of the quaternion.
 * - @quat.y: The Y component of the quaternion.
 * - @quat.z: The Z component of the quaternion.
 * - @quat.w: The W (real) component of the quaternion.
 *
 * The resulting matrix applies the rotation described by the quaternion
 * when multiplied by a vector or another matrix.
 *
 * Note: The function assumes the matrix `mat4` is a 4x4 matrix, and it resets
 *       the translation components to zero, leaving the rotation part intact.
 */
void CMatrix4Df::InitRotateTransform(const glm::quat& sQuat)
{
	glm::mat4 mat = glm::mat4_cast(sQuat);

	CMatrix4Df res(mat);

	*this = res;
}

void CMatrix4Df::InitRotateTransformDir(const SVector3Df& vDir)
{
	SVector3Df vUp(0.0f, 1.0f, 0.0f);
	InitCameraTransform(vDir, vUp);
}

void CMatrix4Df::InitRotateTransformDir(const glm::vec3& vDir)
{
	glm::vec3 vUp(0.0f, 1.0f, 0.0f);
	InitCameraTransform(vDir, vUp);
}
