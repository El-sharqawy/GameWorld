#pragma once

#include "vectors.h"

/**
 * Structure representing a quaternion for 3D rotations
 *
 * This structure is used to represent a quaternion, which is commonly used in 3D graphics
 * and physics to perform rotations. Quaternions are a more efficient and stable alternative
 * to Euler angles and matrices for rotating objects in 3D space, avoiding issues such as gimbal lock.
 *
 * @x: The x component of the quaternion.
 * @y: The y component of the quaternion.
 * @z: The z component of the quaternion.
 * @w: The w component of the quaternion (scalar part).
 */
struct SQuaternion
{
	/* The components of the quaternion: x, y, z (vector part), and w (scalar part) */
	float x, y, z, w;

	/**
	 * Default constructor, initializes all components to zero.
	 */
	SQuaternion() = default;

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
	SQuaternion(const float fAngle, const SVector3Df& vec);

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
	SQuaternion(const float _x, const float _y, const float _z, const float _w);

	/**
	 * Computes the length (magnitude) of the quaternion.
	 *
	 * This function calculates the length of the quaternion using the formula:
	 * length = sqrt(x^2 + y^2 + z^2 + w^2).
	 *
	 * @return The magnitude of the quaternion.
	 */
	float length() const;

	/**
	 * Normalizes the quaternion to unit length.
	 *
	 * This function normalizes the quaternion by dividing each of its components by its length.
	 * Normalization ensures that the quaternion represents a unit vector, which is necessary
	 * for many operations in 3D rotations.
	 */
	void normalize();

	/**
	 * SQuaternion::conjugate() const - Returns the conjugate of the quaternion.
	 *
	 * The conjugate of a quaternion is the result of negating the vector part (x, y, z) while keeping
	 * the scalar part (w) unchanged. This operation is useful for finding the inverse of a quaternion.
	 *
	 * @return A new quaternion representing the conjugate of the current quaternion.
	 */
	SQuaternion conjugate() const;

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
	SVector3Df ToDegrees(bool bSafety = false) const;

	/**
	 * Checks if the quaternion represents no rotation.
	 *
	 * This function checks if the quaternion is a zero quaternion, meaning it represents no rotation
	 * (i.e., all components are zero). A zero quaternion does not perform any rotation and is used
	 * for initialization or invalid states.
	 *
	 * @return true if the quaternion is zero (no rotation), false otherwise.
	 */
	bool IsZero() const;
};

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
SQuaternion operator*(const SQuaternion& QuatLeft, const SQuaternion& QuatRight);

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
SQuaternion operator*(const SQuaternion& Quat, const SVector3Df& vec);
