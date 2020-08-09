#include "stdafx.h"

//given a matrix, it will reduce the values in all rows under the matrix[rowIndex][rowIndex] value to 0. 
//so the column below matrix[rowIndex][rowIndex] is made 0.
//if row reduction is done for all rows, it makes the lower diagonal entries of the matrix 0
void rowReduction(double **matrix, int rowIndex, int rows, int columns) {
	int numberOfRowsToReduce = rows - 1 - rowIndex;

	double rowIndexDataValue = matrix[rowIndex][rowIndex];

	for (size_t i = 1; i <= numberOfRowsToReduce; i++) {
		double rowValueToReduce = matrix[rowIndex + i][rowIndex];
		double multiplicant = 0;
		if (rowIndexDataValue != 0)
			multiplicant = rowValueToReduce / rowIndexDataValue;

		//multiply the multiplicant rowindex row values and substract it from row to be reduced
		for (size_t j = 0; j < columns; j++) {
			matrix[rowIndex + i][j] = matrix[rowIndex + i][j] - multiplicant * matrix[rowIndex][j];
		}
	}
	////debug
	//printf("\nRow reduced matrix is:\n");
	//printMatrix(matrix, rows, columns);
	//printf("\n");
}


//a method of calculating determinants. Much faster than naive determinant evaluation but works only for square matrices
double PivotalCondensationMethod(double **matrix, int rows, int columns) {
	if (rows != columns) {
		printf("Only square matrices can evaluated using pivotal condenstaion method");
	}
	if (rows == 1 && columns == 1) {
		return matrix[0][0];
	}
	double v = matrix[0][0];
	rowReduction(matrix, 0, rows, columns);
	double **reducedMatrix = viewOfMatrix(matrix, rows, columns, 0);

	////////debug
	//printf("%lf * \n", v);
	//printMatrix(reducedMatrix, rows - 1, columns - 1);
	//printf("\n");

	double pivotalValue = v* PivotalCondensationMethod(reducedMatrix, rows - 1, columns - 1);
	//deallocate the reduced matrix for matrix[0][col] element as we wont be using the same matrix again
	for (int i = 0; i < rows - 1; i++)
		free(reducedMatrix[i]);
	free(reducedMatrix);
	return pivotalValue;
}


//it elimnates the first row of the matrix and a given column col from the matrix to give a filtered view
//it is the responsibility of the user of this method to free the memory used by the view
//todo: parameterize the row to be deleted as well like col.
double** viewOfMatrix(double **matrix, int rows, int columns, int col) {
	double **reducedMatrix;
	reducedMatrix = (double **)malloc((rows - 1) * sizeof(double *));
	for (int i = 0; i < rows - 1; i++) {
		reducedMatrix[i] = (double*)malloc((columns - 1) * sizeof(double));
	}

	//skip the first row which is the column of the matrix[0][col] element
	for (size_t i = 1; i < rows; i++) {
		int colIndex = 0;//the column value for the reduced view matrix would be different from the actual matrix. So define your own col index variable
		for (size_t j = 0; j < columns; j++) {
			//and skip the column of the  matrix[0][col] element
			if (j != col) {
				reducedMatrix[i - 1][colIndex] = matrix[i][j];
				colIndex++;
			}
		}
	}
	return reducedMatrix;
}

//naive determinant implementation that is slow compared to eigen library.
//at each step, the matrix will remain square
double Determinant(double **matrix, int rows, int columns) {
	if (rows != columns) {
		printf("Determinant can only be calculated for square matices");
		return NAN;
	}

	if (rows == 1 && columns == 1) {
		return matrix[0][0];
	}

	if (rows == 2 && columns == 2) {
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	}

	int sign = -1;
	double det = 0;
	for (size_t col = 0; col < columns; col++)
	{
		sign = sign * -1;
		//we have to create a view of the matrix excluding the row and column of the matrix[0][col] element
		//now the row to exclude would always be the first row of the view
		double **reducedMatrix = viewOfMatrix(matrix, rows, columns, col);


		////for debug
		/*printf("\n");
		printf("%lf * %d * \n", matrix[0][col], sign);
		printMatrix(reducedMatrix, rows - 1, columns - 1);*/

		det = det + matrix[0][col] * Determinant(&reducedMatrix[0], rows - 1, columns - 1) * sign;

		//deallocate the reduced matrix for matrix[0][col] element as we wont be using the same matrix again
		for (int i = 0; i < rows - 1; i++)
			free(reducedMatrix[i]);
		free(reducedMatrix);
	}
	return det;
}

//double Determinant(double **matrix, int rows, int columns) {
//	if (rows != columns) {
//		printf("Determinant can only be calculated for square matices");
//		return NAN;
//	}
//
//	if (rows == 1 && columns == 1) {
//		return matrix[0][0];
//	}
//
//	double det = calcDeterminant(matrix, rows, columns);
//	printf("The determinant value of the matrix is: %lf\n", det);
//	return det;
//}

