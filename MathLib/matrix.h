#pragma once

#include <glm/glm.hpp>
#include "utils.h"
#include "vectors.h"
#include "quaternion.h"

/**
 * CMatrix4Df: A 4x4 matrix class.
 *
 * This class represents a 4x4 matrix, commonly used in linear algebra and 3D graphics.
 * It provides various constructors for initialization and supports common matrix operations.
 *
 * Members:
 *   - mat4: A 4x4 array of floating-point values representing the matrix elements.
 */
class CMatrix4Df
{
public:
	float mat4[4][4];

	/**
	 * Default constructor, initializes all elements to zero.
	 */
	CMatrix4Df() = default;

	/**
	 * Constructor that initializes the matrix with specific values.
	 *
	 * @param a00, a01, ..., a33: Individual matrix elements.
	 */
	CMatrix4Df(const float a00, const float a01, const float a02, const float a03,
		const float a10, const float a11, const float a12, const float a13,
		const float a20, const float a21, const float a22, const float a23,
		const float a30, const float a31, const float a32, const float a33)
	{
		mat4[0][0] = a00; mat4[0][1] = a01; mat4[0][2] = a02; mat4[0][3] = a03;
		mat4[1][0] = a10; mat4[1][1] = a11; mat4[1][2] = a12; mat4[1][3] = a13;
		mat4[2][0] = a20; mat4[2][1] = a21; mat4[2][2] = a22; mat4[2][3] = a23;
		mat4[3][0] = a20; mat4[3][1] = a31; mat4[3][2] = a32; mat4[3][3] = a33;
	}
	
	/**
	 * Constructor that initializes the matrix from a GLM matrix.
	 *
	 * @param glmMat: A GLM matrix to copy from.
	 */
	CMatrix4Df(const glm::mat4& glmMat)
	{
		mat4[0][0] = glmMat[0][0]; mat4[0][1] = glmMat[0][1]; mat4[0][2] = glmMat[0][2]; mat4[0][3] = glmMat[0][3];
		mat4[1][0] = glmMat[1][0]; mat4[1][1] = glmMat[1][1]; mat4[1][2] = glmMat[1][2]; mat4[1][3] = glmMat[1][3];
		mat4[2][0] = glmMat[2][0]; mat4[2][1] = glmMat[2][1]; mat4[2][2] = glmMat[2][2]; mat4[2][3] = glmMat[2][3];
		mat4[3][0] = glmMat[3][0]; mat4[3][1] = glmMat[3][1]; mat4[3][2] = glmMat[3][2]; mat4[3][3] = glmMat[3][3];
	}

	/**
	 * Copy constructor.
	 *
	 * @param cMat: The matrix to copy.
	 */
	CMatrix4Df(const CMatrix4Df& cMat)
	{
		mat4[0][0] = cMat.mat4[0][0]; mat4[0][1] = cMat.mat4[0][1]; mat4[0][2] = cMat.mat4[0][2]; mat4[0][3] = cMat.mat4[0][3];
		mat4[1][0] = cMat.mat4[1][0]; mat4[1][1] = cMat.mat4[1][1]; mat4[1][2] = cMat.mat4[1][2]; mat4[1][3] = cMat.mat4[1][3];
		mat4[2][0] = cMat.mat4[2][0]; mat4[2][1] = cMat.mat4[2][1]; mat4[2][2] = cMat.mat4[2][2]; mat4[2][3] = cMat.mat4[2][3];
		mat4[3][0] = cMat.mat4[3][0]; mat4[3][1] = cMat.mat4[3][1]; mat4[3][2] = cMat.mat4[3][2]; mat4[3][3] = cMat.mat4[3][3];
	}

	/**
	 * Copy constructor.
	 *
	 * @param vec1: A 4D Vector to copy from first row.
	 */
	CMatrix4Df(const SVector4Df& vec1, const SVector4Df& vec2, const SVector4Df& vec3, const SVector4Df& vec4)
	{
		mat4[0][0] = vec1.x; mat4[0][1] = vec1.y; mat4[0][2] = vec1.z; mat4[0][3] = vec1.w;
		mat4[1][0] = vec2.x; mat4[1][1] = vec2.y; mat4[1][2] = vec2.z; mat4[1][3] = vec2.w;
		mat4[2][0] = vec3.x; mat4[2][1] = vec3.y; mat4[2][2] = vec3.z; mat4[2][3] = vec3.w;
		mat4[3][0] = vec4.x; mat4[3][1] = vec4.y; mat4[3][2] = vec4.z; mat4[3][3] = vec4.w;
	}

	/**
	 * Multiplies two CMatrix4Df matrices.
	 *
	 * This operator overloads the multiplication operator for CMatrix4Df objects,
	 * performing matrix multiplication.
	 *
	 * @param rightMat: The right-hand side matrix.
	 *
	 * @return The product of the two matrices.
	 */
	CMatrix4Df operator*(const CMatrix4Df& rightMat)
	{
		CMatrix4Df newMat{};
		for (int8_t i = 0; i < 4; i++)
		{
			for (int8_t j = 0; j < 4; j++)
			{
				newMat.mat4[i][j] = mat4[i][0] * rightMat.mat4[0][j] + mat4[i][1] * rightMat.mat4[1][j] + mat4[i][2] * rightMat.mat4[2][j] + mat4[i][3] * rightMat.mat4[3][j];
			}
		}
		return (newMat);

	}

	/**
	 * Multiplies a CMatrix4Df and an SVector4Df.
	 *
	 * This operator overloads the multiplication operator for CMatrix4Df and SVector4Df objects,
	 * performing matrix-vector multiplication.
	 *
	 * @param vec: The vector to multiply.
	 *
	 * @return The product of the matrix and the vector.
	 */
	SVector4Df operator*(const SVector4Df& vec)
	{
		SVector4Df newVec{};

		newVec.x = mat4[0][0] * vec.x + mat4[0][1] * vec.y + mat4[0][2] * vec.z + mat4[0][3] * vec.w;
		newVec.y = mat4[1][0] * vec.x + mat4[1][1] * vec.y + mat4[1][2] * vec.z + mat4[1][3] * vec.w;
		newVec.z = mat4[2][0] * vec.x + mat4[2][1] * vec.y + mat4[2][2] * vec.z + mat4[2][3] * vec.w;
		newVec.w = mat4[3][0] * vec.x + mat4[3][1] * vec.y + mat4[3][2] * vec.z + mat4[3][3] * vec.w;
	}

	CMatrix4Df operator*(float scalar)
	{
		CMatrix4Df result;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				result.mat4[i][j] = mat4[i][j] * scalar;
			}
		}
		return result;
	}

	/**
	 * Provides a pointer to the underlying matrix data.
	 *
	 * This operator overload allows direct access to the matrix elements as a float array.
	 *
	 * @return A pointer to the first element of the matrix.
	 */
	operator const float* () const
	{
		return (&(mat4[0][0]));
	}

	/**
	 * Prints the matrix elements to the console.
	 *
	 * This function prints the matrix elements in a readable format.
	 */
	void print() const
	{
		for (int8_t i = 0; i < 4; i++)
		{
			printf("%f - %f - %f - %f", mat4[i][0], mat4[i][1], mat4[i][2], mat4[i][3]);
		}
	}

	/**
	 * Accesses the underlying data of the matrix.
	 *
	 * - This function provides a pointer to the first element of the
	 * private `mat4` 4x4 matrix, allowing direct manipulation of the matrix data.
	 *
	 * @return A pointer to the first element of the `mat4` matrix.
	 */
	float* data()
	{
		return (&(mat4[0][0]));
	}

	/**
	 * returns a reference to the private `mat4` matrix, It allows read-only access to the 4x4 matrix stored in the class
	 *
	 * @return A constant reference to the `mat4` matrix.
	 */
	const float(&GetMatrix() const)[4][4]
	{
		return mat4;
	}

	/**
	 * Takes a 4x4 array of floats and updates the private `mat4` matrix with the provided values.
	 *
	 * @param values: A 4x4 array containing the new values for the matrix.
	 */
	void SetMatrix(const float values[4][4])
	{
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				mat4[i][j] = values[i][j];
			}
		}
	}

	/**
	 * Transposes the matrix.
	 *
	 * This function returns a new matrix that is the transpose of the current matrix.
	 *
	 * @return The transposed matrix.
	 */
	CMatrix4Df Transpose() const
	{
		CMatrix4Df newMat{};

		for (int8_t i = 0; i < 4; i++)
		{
			for (int8_t j = 0; j < 4; j++)
			{
				newMat.mat4[i][j] = mat4[j][i];
			}
		}

		return (newMat);
	}

	/**
	 * Calculates the determinant of the matrix.
	 *
	 * This function calculates the determinant of the 4x4 matrix using a direct formula.
	 *
	 * @return The determinant of the matrix.
	 */
	float Determinant() const
	{
		float fDet = 
			  mat4[0][0] * mat4[1][1] * mat4[2][2] * mat4[3][3] - mat4[0][0] * mat4[1][1] * mat4[2][3] * mat4[3][2] + mat4[0][0] * mat4[1][2] * mat4[2][3] * mat4[3][1] - mat4[0][0] * mat4[1][2] * mat4[2][1] * mat4[3][3]
			+ mat4[0][0] * mat4[1][3] * mat4[2][1] * mat4[3][2] - mat4[0][0] * mat4[1][3] * mat4[2][2] * mat4[3][1] - mat4[0][1] * mat4[1][2] * mat4[2][3] * mat4[3][0] + mat4[0][1] * mat4[1][2] * mat4[2][0] * mat4[3][3]
			- mat4[0][1] * mat4[1][3] * mat4[2][0] * mat4[3][2] + mat4[0][1] * mat4[1][3] * mat4[2][2] * mat4[3][0] - mat4[0][1] * mat4[1][0] * mat4[2][2] * mat4[3][3] + mat4[0][1] * mat4[1][0] * mat4[2][3] * mat4[3][2]
			+ mat4[0][2] * mat4[1][3] * mat4[2][0] * mat4[3][1] - mat4[0][2] * mat4[1][3] * mat4[2][1] * mat4[3][0] + mat4[0][2] * mat4[1][0] * mat4[2][1] * mat4[3][3] - mat4[0][2] * mat4[1][0] * mat4[2][3] * mat4[3][1]
			+ mat4[0][2] * mat4[1][1] * mat4[2][3] * mat4[3][0] - mat4[0][2] * mat4[1][1] * mat4[2][0] * mat4[3][3] - mat4[0][3] * mat4[1][0] * mat4[2][1] * mat4[3][2] + mat4[0][3] * mat4[1][0] * mat4[2][2] * mat4[3][1]
			- mat4[0][3] * mat4[1][1] * mat4[2][2] * mat4[3][0] + mat4[0][3] * mat4[1][1] * mat4[2][0] * mat4[3][2] - mat4[0][3] * mat4[1][2] * mat4[2][0] * mat4[3][2] + mat4[0][3] * mat4[1][2] * mat4[2][1] * mat4[3][0];

		return (fDet);
	}

	/**
	 * Calculates the determinant of a 3x3 submatrix, Calcualted as Same as GLM Function.
	 *
	 * This function calculates the determinant of a 3x3 submatrix, which is a crucial
	 * step in the cofactor expansion method for computing the determinant of a 4x4 matrix.
	 *
	 * @return The determinant of the 3x3 submatrix.
	 */
	float DeterminantSub() const
	{
		float SubFactor00 = mat4[2][2] * mat4[3][3] - mat4[3][2] * mat4[2][3];
		float SubFactor01 = mat4[2][1] * mat4[3][3] - mat4[3][1] * mat4[2][3];
		float SubFactor02 = mat4[2][1] * mat4[3][2] - mat4[3][1] * mat4[2][2];
		float SubFactor03 = mat4[2][0] * mat4[3][3] - mat4[3][0] * mat4[2][3];
		float SubFactor04 = mat4[2][0] * mat4[3][2] - mat4[3][0] * mat4[2][2];
		float SubFactor05 = mat4[2][0] * mat4[3][1] - mat4[3][0] * mat4[2][1];

		SVector4Df DetCof(
			+ (mat4[1][1] * SubFactor00 - mat4[1][2] * SubFactor01 + mat4[1][3] * SubFactor02),
			- (mat4[1][0] * SubFactor00 - mat4[1][2] * SubFactor03 + mat4[1][3] * SubFactor04),
			+ (mat4[1][0] * SubFactor01 - mat4[1][1] * SubFactor03 + mat4[1][3] * SubFactor05),
			- (mat4[1][0] * SubFactor02 - mat4[1][1] * SubFactor04 + mat4[1][2] * SubFactor05)
		);

		return mat4[0][0] * DetCof[0] + mat4[0][1] * DetCof[1] + mat4[0][2] * DetCof[2] + mat4[0][3] * DetCof[3];
	}

	/**
	 * Calculates the inverse of the matrix.
	 *
	 * This function calculates the inverse of the 4x4 matrix using the adjugate matrix and the determinant.
	 *
	 * @return The inverse of the matrix.
	 */
	CMatrix4Df Inverse() const
	{
		// Compute the reciprocal determinant first

		float det = Determinant();

		if (det == 0.0f)
		{
			ASSERT(det == 0.0f, "Matrix Determinant Is 0");
			return (*this);
		}

		// Calculate the inverse determinant
		float fInvDet = 1.0f / det;

		// Calculate the adjugate matrix (transpose of the cofactor matrix) (https://byjus.com/maths/inverse-matrix/)
		CMatrix4Df res{};

		// Calcualte First Element --- Inverse Is Positive
		res.mat4[0][0] = fInvDet * (mat4[1][1] * (mat4[2][2] * mat4[3][3] - mat4[2][3] * mat4[3][2]) + mat4[1][2] * 
			(mat4[2][3] * mat4[3][1] - mat4[2][1] * mat4[3][3]) + mat4[1][3] * (mat4[2][1] * mat4[3][2] - mat4[2][2] * mat4[3][1]));

		// Same as the previous just switch the [1][x] to [0][x] --- Inverse Is Negative
		res.mat4[0][1] = -fInvDet * (mat4[0][1] * (mat4[2][2] * mat4[3][3] - mat4[2][3] * mat4[3][2]) + mat4[0][2] *
			(mat4[2][3] * mat4[3][1] - mat4[2][1] * mat4[3][3]) + mat4[0][3] * (mat4[2][1] * mat4[3][2] - mat4[2][2] * mat4[3][1]));

		// Same as the previous just switch the [2][x] to [1][x] --- Inverse Is Positive
		res.mat4[0][2] = fInvDet * (mat4[0][1] * (mat4[1][2] * mat4[3][3] - mat4[1][3] * mat4[3][2]) + mat4[0][2] *
			(mat4[1][3] * mat4[3][1] - mat4[1][1] * mat4[3][3]) + mat4[0][3] * (mat4[1][1] * mat4[3][2] - mat4[1][2] * mat4[3][1]));

		// Same as the previous just switch the [3][x] to [2][x] --- Inverse Is Negative
		res.mat4[0][3] = -fInvDet * (mat4[0][1] * (mat4[1][2] * mat4[2][3] - mat4[1][3] * mat4[3][2]) + mat4[0][2] *
			(mat4[1][3] * mat4[2][1] - mat4[1][1] * mat4[2][3]) + mat4[0][3] * (mat4[1][1] * mat4[2][2] - mat4[1][2] * mat4[2][1]));

		// Same as the FIRST just switch the [x][1] to [x][0] --- Inverse Is Negative
		res.mat4[1][0] = -fInvDet * (mat4[1][0] * (mat4[2][2] * mat4[3][3] - mat4[2][3] * mat4[3][2]) + mat4[1][2] *
			(mat4[2][3] * mat4[3][0] - mat4[2][0] * mat4[3][3]) + mat4[1][3] * (mat4[2][0] * mat4[3][2] - mat4[2][2] * mat4[3][0]));

		// Same as the previous just switch the [1][x] to [0][x] --- Inverse Is Positive
		res.mat4[1][1] = fInvDet * (mat4[0][0] * (mat4[2][2] * mat4[3][3] - mat4[2][3] * mat4[3][2]) + mat4[0][2] *
			(mat4[2][3] * mat4[3][0] - mat4[2][0] * mat4[3][3]) + mat4[0][3] * (mat4[2][0] * mat4[3][2] - mat4[2][2] * mat4[3][0]));

		// Same as the previous just switch the [2][x] to [1][x] --- Inverse Is Negative
		res.mat4[1][2] = -fInvDet * (mat4[0][0] * (mat4[1][2] * mat4[3][3] - mat4[1][3] * mat4[3][2]) + mat4[0][2] *
			(mat4[1][3] * mat4[3][0] - mat4[1][0] * mat4[3][3]) + mat4[0][3] * (mat4[1][0] * mat4[3][2] - mat4[1][2] * mat4[3][0]));

		// Same as the previous just switch the [3][x] to [2][x] --- Inverse Is Positive
		res.mat4[1][3] = fInvDet * (mat4[0][0] * (mat4[1][2] * mat4[2][3] - mat4[1][3] * mat4[2][2]) + mat4[0][2] *
			(mat4[1][3] * mat4[2][0] - mat4[1][0] * mat4[2][3]) + mat4[0][3] * (mat4[1][0] * mat4[2][2] - mat4[1][2] * mat4[2][0]));

		// Same as [1][0] ELEMENT calculation, just switch the [x][2] to [x][1] --- Inverse Is Positive
		res.mat4[2][0] = fInvDet * (mat4[1][0] * (mat4[2][1] * mat4[3][3] - mat4[2][3] * mat4[3][1]) + mat4[1][1] *
			(mat4[2][3] * mat4[3][0] - mat4[2][0] * mat4[3][3]) + mat4[1][3] * (mat4[2][0] * mat4[3][1] - mat4[2][1] * mat4[3][0]));

		// Same as the previous just switch the [1][x] to [0][x] --- Inverse Is Negative
		res.mat4[2][1] = -fInvDet * (mat4[0][0] * (mat4[2][1] * mat4[3][3] - mat4[2][3] * mat4[3][1]) + mat4[0][1] *
			(mat4[2][3] * mat4[3][0] - mat4[2][0] * mat4[3][3]) + mat4[0][3] * (mat4[2][0] * mat4[3][1] - mat4[2][1] * mat4[3][0]));

		// Same as the previous just switch the [2][x] to [1][x] --- Inverse Is Positive
		res.mat4[2][2] = fInvDet * (mat4[0][0] * (mat4[1][1] * mat4[3][3] - mat4[1][3] * mat4[3][1]) + mat4[0][1] *
			(mat4[1][3] * mat4[3][0] - mat4[1][0] * mat4[3][3]) + mat4[0][3] * (mat4[1][0] * mat4[3][1] - mat4[1][1] * mat4[3][0]));

		// Same as the previous just switch the [3][x] to [2][x] --- Inverse Is Negative
		res.mat4[2][3] = -fInvDet * (mat4[0][0] * (mat4[1][1] * mat4[2][3] - mat4[1][3] * mat4[2][1]) + mat4[0][1] *
			(mat4[1][3] * mat4[2][0] - mat4[1][0] * mat4[2][3]) + mat4[0][3] * (mat4[1][0] * mat4[2][1] - mat4[1][1] * mat4[2][0]));

		// Same as [2][0] ELEMENT calculation, just switch the [x][3] to [x][2] --- Inverse Is Negative
		res.mat4[3][0] = -fInvDet * (mat4[1][0] * (mat4[2][1] * mat4[3][2] - mat4[2][2] * mat4[3][1]) + mat4[1][1] *
			(mat4[2][2] * mat4[3][0] - mat4[2][0] * mat4[3][2]) + mat4[1][2] * (mat4[2][0] * mat4[3][1] - mat4[2][1] * mat4[3][0]));

		// Same as the previous just switch the [1][x] to [0][x] --- Inverse Is Positive
		res.mat4[3][1] = fInvDet * (mat4[0][0] * (mat4[2][1] * mat4[3][2] - mat4[2][2] * mat4[3][1]) + mat4[0][1] *
			(mat4[2][2] * mat4[3][0] - mat4[2][0] * mat4[3][2]) + mat4[0][2] * (mat4[2][0] * mat4[3][1] - mat4[2][1] * mat4[3][0]));

		// Same as the previous just switch the [2][x] to [1][x] --- Inverse Is Negative
		res.mat4[3][2] = -fInvDet * (mat4[0][0] * (mat4[1][1] * mat4[3][2] - mat4[1][2] * mat4[3][1]) + mat4[0][1] *
			(mat4[1][2] * mat4[3][0] - mat4[1][0] * mat4[3][2]) + mat4[0][2] * (mat4[2][0] * mat4[3][1] - mat4[1][1] * mat4[3][0]));

		// Same as the previous just switch the [3][x] to [2][x] --- Inverse Is Positive
		res.mat4[3][3] = fInvDet * (mat4[0][0] * (mat4[1][1] * mat4[2][2] - mat4[1][2] * mat4[2][1]) + mat4[0][1] *
			(mat4[1][2] * mat4[3][0] - mat4[1][0] * mat4[2][2]) + mat4[0][2] * (mat4[2][0] * mat4[2][1] - mat4[1][1] * mat4[2][0]));

		// The Pattern -> Previous Calculation (make sure of fInvDet sign) and each time make [x][n] = [x-1][n] , ex: [3][2] -> [2][2], [2][2] -> [1][2], [1][2] -> [0][2]
		return (res);
	}

	/**
	 * Calculates the inverse of the matrix, Calcualted as Same as GLM Function.
	 *
	 * This function calculates the inverse of the 4x4 matrix using the adjugate matrix and the determinant.
	 *
	 * @return The inverse of the matrix.
	 */
	CMatrix4Df InverseSub() const
	{
		float Coef00 = mat4[2][2] * mat4[3][3] - mat4[3][2] * mat4[2][3];
		float Coef02 = mat4[1][2] * mat4[3][3] - mat4[3][2] * mat4[1][3];
		float Coef03 = mat4[1][2] * mat4[2][3] - mat4[2][2] * mat4[1][3];

		float Coef04 = mat4[2][1] * mat4[3][3] - mat4[3][1] * mat4[2][3];
		float Coef06 = mat4[1][1] * mat4[3][3] - mat4[3][1] * mat4[1][3];
		float Coef07 = mat4[1][1] * mat4[2][3] - mat4[2][1] * mat4[1][3];

		float Coef08 = mat4[2][1] * mat4[3][2] - mat4[3][1] * mat4[2][3];
		float Coef10 = mat4[1][1] * mat4[3][2] - mat4[3][1] * mat4[1][2];
		float Coef11 = mat4[1][1] * mat4[2][2] - mat4[2][1] * mat4[1][2];

		float Coef12 = mat4[2][0] * mat4[3][3] - mat4[3][0] * mat4[2][3];
		float Coef14 = mat4[1][0] * mat4[3][3] - mat4[3][0] * mat4[1][3];
		float Coef15 = mat4[1][0] * mat4[2][3] - mat4[2][0] * mat4[1][3];

		float Coef16 = mat4[2][0] * mat4[3][2] - mat4[3][0] * mat4[2][2];
		float Coef18 = mat4[1][0] * mat4[3][2] - mat4[3][0] * mat4[1][2];
		float Coef19 = mat4[1][0] * mat4[2][2] - mat4[2][0] * mat4[1][2];

		float Coef20 = mat4[2][0] * mat4[3][1] - mat4[3][0] * mat4[2][1];
		float Coef22 = mat4[1][0] * mat4[3][1] - mat4[3][0] * mat4[1][1];
		float Coef23 = mat4[1][0] * mat4[2][1] - mat4[2][0] * mat4[1][1];

		SVector4Df Fac0(Coef00, Coef00, Coef02, Coef03);
		SVector4Df Fac1(Coef04, Coef04, Coef06, Coef07);
		SVector4Df Fac2(Coef08, Coef08, Coef10, Coef11);
		SVector4Df Fac3(Coef12, Coef12, Coef14, Coef15);
		SVector4Df Fac4(Coef16, Coef16, Coef18, Coef19);
		SVector4Df Fac5(Coef20, Coef20, Coef22, Coef23);


		SVector4Df Vec0(mat4[1][0], mat4[0][0], mat4[0][0], mat4[0][0]);
		SVector4Df Vec1(mat4[1][1], mat4[0][1], mat4[0][1], mat4[0][1]);
		SVector4Df Vec2(mat4[1][2], mat4[0][2], mat4[0][2], mat4[0][2]);
		SVector4Df Vec3(mat4[1][3], mat4[0][3], mat4[0][3], mat4[0][3]);

		SVector4Df Inv0(Vec1 * Fac0 - Vec2 * Fac1 + Vec3 * Fac2);
		SVector4Df Inv1(Vec0 * Fac0 - Vec2 * Fac3 + Vec3 * Fac4);
		SVector4Df Inv2(Vec0 * Fac1 - Vec1 * Fac3 + Vec3 * Fac5);
		SVector4Df Inv3(Vec0 * Fac2 - Vec1 * Fac4 + Vec2 * Fac5);

		SVector4Df SignA(+1, -1, +1, -1);
		SVector4Df SignB(-1, +1, -1, +1);

		CMatrix4Df Inverse(Inv0 * SignA, Inv1 * SignB, Inv2 * SignA, Inv3 * SignB);

		SVector4Df Row0(Inverse.mat4[0][0], Inverse.mat4[1][0], Inverse.mat4[2][0], Inverse.mat4[3][0]);

		SVector4Df Dot0(mat4[0] * Row0);

		float Dot1 = (Dot0.x + Dot0.y) + (Dot0.z + Dot0.w);

		float OneOverDeterminant = 1.0f / Dot1;

		return Inverse * OneOverDeterminant;
	}

	/*
	 * Calcualte Matrix Inverse Using Gauss Jordan Method.
	 */
	CMatrix4Df InverseGJ(const CMatrix4Df& inputMatrix)
	{
		CMatrix4Df inverseMatrix{};
		float augmented[4][8] = { 0 };

		// Get the input matrix
		const auto& mat = inputMatrix.GetMatrix();

		// Create the augmented matrix
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				augmented[i][j] = mat[i][j]; // Copy the original matrix
			}
			augmented[i][i + 4] = 1; // Add the identity matrix
		}

		// Perform Gaussian elimination
		for (int i = 0; i < 4; i++)
		{
			// Normalize the current row
			float diagElement = augmented[i][i];
			for (int j = 0; j < 8; j++)
			{
				augmented[i][j] /= diagElement;
			}

			// Eliminate other rows
			for (int k = 0; k < 4; k++)
			{
				if (k != i)
				{
					float factor = augmented[k][i];
					for (int j = 0; j < 8; j++)
					{
						augmented[k][j] -= factor * augmented[i][j];
					}
				}
			}
		}

		// Extract the inverse matrix
		float result[4][4] = { 0 };
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result[i][j] = augmented[i][j + 4];
			}
		}

		// Set the result to the inverse matrix object
		inverseMatrix.SetMatrix(result);
		return inverseMatrix;
	}


	/**
 * Initializes the matrix to zero using memset.
 *
 * This function efficiently sets all elements of the matrix to zero using the `memset` function.
 */
	void InitMemZero()
	{
		arr_mem_zero_ref(mat4);
	}

	/**
	 * InitZero: Initializes the matrix to zero element-wise.
	 *
	 * This function sets all elements of the matrix to zero using a nested loop.
	 */
	void InitZero()
	{
		for (int8_t i = 0; i < 4; i++)
		{
			for (int8_t j = 0; j < 4; j++)
			{
				mat4[i][j] = 0.0f;
			}
		}
	}

	/**
	 * Initializes the matrix with a specified number.
	 *
	 * This function initializes the matrix elements by the specified `fNum` value.
	 *
	 * @param fNum: The value for each matrix element.
	 */
	void InitNum(const float fNum)
	{
		for (int8_t i = 0; i < 4; i++)
		{
			for (int8_t j = 0; j < 4; j++)
			{
				mat4[i][j] = fNum;
			}
		}
	}

	/**
	 * Initializes the matrix to the identity matrix.
	 *
	 * This function sets the matrix to the identity matrix, which has ones on the
	 * diagonal and zeros elsewhere.
	 */
	void InitIdentity()
	{
		mat4[0][0] = 1.0f; mat4[0][1] = 0.0f; mat4[0][2] = 0.0f; mat4[0][3] = 0.0f;
		mat4[1][0] = 0.0f; mat4[1][1] = 1.0f; mat4[1][2] = 0.0f; mat4[1][3] = 0.0f;
		mat4[2][0] = 0.0f; mat4[2][1] = 0.0f; mat4[2][2] = 1.0f; mat4[2][3] = 0.0f;
		mat4[3][0] = 0.0f; mat4[3][1] = 0.0f; mat4[3][2] = 0.0f; mat4[3][3] = 1.0f;
	}

	void InitScaleTransform(const float fScaleX, const float fScaleY, const float fScaleZ);
	void InitScaleTransform(const float fScale);
	void InitScaleTransform(const SVector3Df& vScale);
	void InitScaleTransform(const glm::vec3& vScale);


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
	void InitRotateTransform(const float fRotateX, const float fRotateY, const float fRotateZ);

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
	void InitRotateTransform(const SVector3Df& vRotation);

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
	void InitRotateTransform(const glm::vec3& vRotation);

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
	void InitRotateTransformZYX(const float fRotateX, const float fRotateY, const float fRotateZ);

	/**
	 * Initializes the matrix for a rotation transformation using a quaternion.
	 * @param sQuat: The quaternion representing the rotation to be applied.
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
	void InitRotateTransform(const SQuaternion& sQuat);

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
	void InitRotateTransform(const glm::quat& sQuat);

	void InitRotateTransformDir(const SVector3Df& vDir);
	void InitRotateTransformDir(const glm::vec3& vDir);

	void InitCameraTransform(const SVector3Df& vTarget, const SVector3Df& vUp);
	void InitCameraTransform(const glm::vec3& vTarget, const glm::vec3& vUp);

protected:
	void InitRoatationX(const float fRotX, const bool bLeftRot = false);
	void InitRoatationY(const float fRotY, const bool bLeftRot = false);
	void InitRoatationZ(const float fRotZ, const bool bLeftRot = false);
};
