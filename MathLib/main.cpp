#include "quaternion.h"
#include "vectors.h"
#include "utils.h"
#include "matrix.h"

int main()
{
	SRANDOM;
	
	glm::mat4 matrix = glm::mat4(
		4.0f, 7.0f, 2.0f, 3.0f,
		3.0f, 6.0f, 1.0f, 4.0f,
		0.0f, 0.0f, 3.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 8.0f
	);

	CMatrix4Df mat4(
		4.0f, 7.0f, 2.0f, 3.0f,
		3.0f, 6.0f, 1.0f, 4.0f,
		0.0f, 0.0f, 3.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 8.0f
	);

	CMatrix4Df inverseMatrix{};
	inverseMatrix = mat4.InverseGJ(mat4);

	printf("_________________________\n");
	for (int i = 0; i < 4; i++)
		printf("%f %f %f %f\n", mat4.mat4[i][0], mat4.mat4[i][1], mat4.mat4[i][2], mat4.mat4[i][3]);

	printf("_________________________\n");
	printf("Determinant: %f\n", mat4.Determinant());
	printf("Determinant: %f\n", mat4.DeterminantSub());
	printf("Determinant: %f\n", glm::determinant(matrix));

	auto matrix2 = glm::inverse(matrix);

	printf("_________________________MY INVERSED MATRIX METHOD_________________________\n");
	for (int i = 0; i < 4; i++)
		printf("%f %f %f %f\n", inverseMatrix.mat4[i][0], inverseMatrix.mat4[i][1], inverseMatrix.mat4[i][2], inverseMatrix.mat4[i][3]);

	printf("_________________________OPEN GL INVERSED MATRIX METHOD_________________________\n");
	for (int i = 0; i < 4; i++)
		printf("%f %f %f %f\n", matrix2[i][0], matrix2[i][1], matrix2[i][2], matrix2[i][3]);

	printf("__________________________________________________\n");
	printf("Determinant: %f\n", inverseMatrix.Determinant());
	printf("Determinant: %f\n", inverseMatrix.DeterminantSub());
	printf("Determinant: %f\n", glm::determinant(matrix2));

	return (EXIT_SUCCESS);
}