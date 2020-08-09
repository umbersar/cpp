#include "stdafx.h"

//mallocs' a 2 dim matrix. The end user of this method should free the memory
double** Create2DimensionalMatrix(int matrixRows, int matrixColumns) {
	double **matA = (double **)malloc(matrixRows * sizeof(double *));
	for (int i = 0; i < matrixRows; i++) {
		matA[i] = (double*)malloc(matrixColumns * sizeof(double));
		//memset(matA[i], 0.0, matrixColumns * sizeof(double));
	}
	return matA;
}

//initializes a 2 dim matrix.
void Init2DimensionalMatrix(double **matrix, int matrixRows, int matrixColumns, double *arr, int sizeOfArray) {
	assert(sizeOfArray == matrixRows*matrixColumns);
	for (size_t i = 0; i < matrixRows; i++) {
		for (size_t j = 0; j < matrixColumns; j++) {
			matrix[i][j] = arr[matrixColumns*i + j];
		}
	}
}

//initializes a 2 dim matrix.
void Init2DimensionalMatrix(double **matrix, int matrixRows, int matrixColumns, double seed) {
	for (size_t i = 0; i < matrixRows; i++) {
		for (size_t j = 0; j < matrixColumns; j++) {
			matrix[i][j] = matrixColumns*i + j + seed;
		}
	}
}

void printMatrix(double **matrix, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			printf("%lf ", matrix[i][j]);
		}
		printf("\n");
	}
}

long factorial(int i) {
	if (i == 0)
		return 1;
	else
		return i * factorial(i - 1);
}

bool sameSign(double x, double y) {
	return ((x < 0) == (y < 0));
}

//this method takes in a function/equation as returns an initial guess of its root.
//taking 0 as the initial guess does not work always. There are some other sohisticated ways to make guesses for the roots of equation.
//this method evaluates f using x and x+1 (or x-1 in case of negative roots) and if the sign of evaluated value of f(x) and f(x+1) is different, then
//the root lies between x and x+1 (or x-1 in case of negative roots) and either value can be used as an initial guesss.
//todo: write it recursively.
//page 299 of C.Xavier book on numerical methods explains this.
double initialGuess(double(*f)(double)) {
	int index = 0;
	double newValue;
	double lastValue = 0;
	double rootGuess = 0;
	while (true) {
		if (index == 500000) {
			//printf("Could not find a reasonable initial guess %d iterations. Relax the max iteration count to see if a solution is possible");
			break;
		}
		newValue = f(rootGuess);
		if (!sameSign(newValue, lastValue))
			return rootGuess;
		rootGuess = rootGuess + 1;
		index++;
	}

	lastValue = 0;
	rootGuess = 0;
	index = 0;
	while (true) {
		if (index == 500000) {
			printf("Could not find a reasonable initial guess. Relax the max iteration count to see if a solution is possible.\n");
			break;
		}
		newValue = f(rootGuess);
		if (!sameSign(newValue, lastValue))
			return rootGuess;
		rootGuess = rootGuess - 1;
		index++;
	}

}