#include "utils.h"
#include "quaternion.h"

/**
 * Constructor that creates a quaternion from an axis of rotation and an angle.
 *
 * This constructor creates a quaternion representing a rotation of the specified angle (in degrees)
 * around the axis defined by the given vector. The angle is divided by 2 and converted to radians
 * to compute the sine and cosine values, which are then used to initialize the quaternion components.
 *
 * @param fAngle: The rotation angle (in degrees).
 * @param vec: The SVector3Df Object vector as axis of rotation.
 */
SQuaternion::SQuaternion(const float fAngle, const SVector3Df& vec)
{
	/* Convert angle to radians and divide by 2 */
	float HalfAngleInRadians = ToRadian(fAngle / 2);

	/* Sine of half the angle */
	float fSinHalfAngle = std::sinf(HalfAngleInRadians);

	/* Cosine of half the angle */
	float fCosHalfAngle = std::cosf(HalfAngleInRadians);

	/* Initialize quaternion components */
	x = vec.x * fSinHalfAngle;
	y = vec.y * fSinHalfAngle;
	z = vec.z * fSinHalfAngle;
	w = fCosHalfAngle;
}

/**
 * Constructor to initialize a quaternion with specific x, y, z, and w values.
 *
 * This constructor initializes the quaternion directly with the provided x, y, z, and w components.
 * It is used when you already have quaternion values to create an instance.
 *
 * @param x: The x component of the quaternion.
 * @param y: The y component of the quaternion.
 * @param z: The z component of the quaternion.
 * @param w: The w (scalar) component of the quaternion.
 */
SQuaternion::SQuaternion(const float _x, const float _y, const float _z, const float _w)
{
	x = _x; /* Set quaternion x component */
	y = _y; /* Set quaternion y component */
	z = _y; /* Set quaternion z component */
	w = _w; /* Set quaternion w (scalar) component */
}

/**
 * Computes the length (magnitude) of the quaternion.
 *
 * This function calculates the length of the quaternion using the formula:
 * length = sqrt(x^2 + y^2 + z^2 + w^2).
 *
 * @return The magnitude of the quaternion.
 */
float SQuaternion::length() const
{
	/* Calculate quaternion magnitude */
	float fLen = sqrtf(x * x + y * y + z * z + w * w);

	/* Return the computed length */
	return (fLen);
}

/**
 * Normalizes the quaternion to unit length.
 *
 * This function normalizes the quaternion by dividing each of its components by its length.
 * Normalization ensures that the quaternion represents a unit vector, which is necessary
 * for many operations in 3D rotations.
 */
void SQuaternion::normalize()
{
	/* Get the magnitude of the quaternion */
	float fLength = length();

	/* Ensure quaternion is not a zero quaternion (avoid division by zero) */
	assert(fLength != 0);

	x /= fLength; /* Normalize x component */
	y /= fLength; /* Normalize y component */
	z /= fLength; /* Normalize z component */
	w /= fLength; /* Normalize w component */
}

/**
 * SQuaternion::conjugate() const - Returns the conjugate of the quaternion.
 *
 * The conjugate of a quaternion is the result of negating the vector part (x, y, z) while keeping
 * the scalar part (w) unchanged. This operation is useful for finding the inverse of a quaternion.
 *
 * @return A new quaternion representing the conjugate of the current quaternion.
 */
SQuaternion SQuaternion::conjugate() const
{
	/* Negate the vector part, keep scalar part unchanged */
	SQuaternion result(-x, -y, -z, w);

	/* Return the conjugate quaternion */
	return (result);
}

/**
 * Converts the quaternion to Euler angles (in degrees).
 *
 * This function converts the quaternion to three Euler angles (pitch, yaw, roll) that represent
 * the same rotation in 3D space. The angles are computed using trigonometric functions.
 * The optional 'bSafety' flag ensures that the conversion handles potential edge cases more safely.
 *
 * @param bSafety: A flag indicating whether to use safe indexing for conversion. Defaults to false.
 * @return A 3D vector containing the Euler angles in degrees (pitch, yaw, roll).
 */
SVector3Df SQuaternion::ToDegrees(bool bSafety) const
{
	/* Initialize an array to store the Euler angles */
	std::array<float, 3> arr{};

	if (!bSafety)
	{
		arr[0] = std::atan2(x * z + y * w, x * w - y * z);	/* Calculate pitch */
		arr[1] = std::acos(-x * x - y * y - z * z - w * w);	/* Calculate yaw */
		arr[2] = std::atan2(x * z - y * w, x * w + y * z);	/* Calculate roll */

		arr[0] = ToDegree(arr[0]); /* Convert to degrees */
		arr[1] = ToDegree(arr[1]); /* Convert to degrees */
		arr[2] = ToDegree(arr[2]); /* Convert to degrees */
	}
	else
	{
		arr.at(0) = std::atan2(x * z + y * w, x * w - y * z);	/* Calculate pitch */
		arr.at(1) = std::acos(-x * x - y * y - z * z - w * w);	/* Calculate yaw */
		arr.at(2) = std::atan2(x * z - y * w, x * w + y * z);	/* Calculate roll */

		arr.at(0) = ToDegree(arr.at(0)); /* Convert to degrees */
		arr.at(1) = ToDegree(arr.at(1)); /* Convert to degrees */
		arr.at(2) = ToDegree(arr.at(2)); /* Convert to degrees */
	}

	/* Return the Euler angles as a vector */
	return (SVector3Df(arr));
}

/**
 * Checks if the quaternion represents no rotation.
 *
 * This function checks if the quaternion is a zero quaternion, meaning it represents no rotation
 * (i.e., all components are zero). A zero quaternion does not perform any rotation and is used
 * for initialization or invalid states.
 *
 * @return true if the quaternion is zero (no rotation), false otherwise.
 */
bool SQuaternion::IsZero() const
{
	return ((x == 0.0f) && (y == 0.0f) && (z == 0.0f) && (w == 0.0f));
}

/**
 * Quaternion multiplication operator.
 *
 * This operator performs multiplication between two quaternions, resulting in a new quaternion that
 * represents the combined rotation of the two input quaternions. Quaternion multiplication is not
 * commutative, meaning the order of multiplication matters.
 *
 * @param QuatLeft: The first quaternion (on the left side of the multiplication).
 * @param QuatRight: The second quaternion (on the right side of the multiplication).
 * @return The result of multiplying the two quaternions.
 */
SQuaternion operator*(const SQuaternion& QuatLeft, const SQuaternion& QuatRight)
{
	float w = (QuatLeft.w * QuatRight.w) - (QuatLeft.x * QuatRight.x) - (QuatLeft.y * QuatRight.y) - (QuatLeft.z * QuatRight.z); /* Compute scalar component */
	float x = (QuatLeft.x * QuatRight.w) + (QuatLeft.w * QuatRight.x) + (QuatLeft.y * QuatRight.z) - (QuatLeft.z * QuatRight.y); /* Compute x component */
	float y = (QuatLeft.y * QuatRight.w) + (QuatLeft.w * QuatRight.y) + (QuatLeft.z * QuatRight.x) - (QuatLeft.x * QuatRight.z); /* Compute y component */
	float z = (QuatLeft.z * QuatRight.w) + (QuatLeft.w * QuatRight.z) + (QuatLeft.x * QuatRight.y) - (QuatLeft.y * QuatRight.x); /* Compute z component */

	/* Create a quaternion with the computed components */
	SQuaternion result(x, y, z, w);

	/* Return the resulting quaternion */
	return (result);
}

/**
 * Quaternion-vector multiplication operator.
 *
 * This operator multiplies a quaternion by a 3D vector. The result is a new vector that represents
 * the original vector rotated by the quaternion. This operation is commonly used to rotate points
 * or directions in 3D space.
 *
 * @param Quat: The quaternion representing the rotation.
 * @param vec: The SVector3Df object to be rotated by the quaternion.
 * @return A new vector that is the result of the rotation.
 */
SQuaternion operator*(const SQuaternion& Quat, const SVector3Df& vec)
{
	float w = (Quat.x * vec.x) - (Quat.y * vec.y) - (Quat.z - vec.z); /* Compute scalar component */
	float x = (Quat.w * vec.x) + (Quat.y * vec.z) - (Quat.z * vec.y); /* Compute x component */
	float y = (Quat.w * vec.y) + (Quat.z * vec.x) - (Quat.x * vec.z); /* Compute y component */
	float z = (Quat.w * vec.z) + (Quat.x * vec.y) - (Quat.y * vec.x); /* Compute z component */

	/* Create a quaternion with the computed components */
	SQuaternion result(x, y, z, w);

	/* Return the resulting quaternion */
	return (result);
}
