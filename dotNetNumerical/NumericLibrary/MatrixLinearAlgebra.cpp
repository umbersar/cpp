#include "stdafx.h"
#include "NumericLibrary.h"

//If A is any square matrix, we can find a constant λ and a vector X such that AX=λX.
//a particular value for λ is called the eigen value and the corresponding vector X is called eigen vector.
//AX=λX can be re-written as AX=λIX. So AX-λIX=0. So [A-λI]X=0 which gives us a coeffecient matrix which is equal to
//subtracting λ from main diagonal entries of A.
//now a non-zero solution X exists for [A-λI]X=0 if and only if the determinant |A-λI|=0
//the solution will give a polynomial equation which can then be solved/roots found using newton-raphson iterative solver
//but before we reach the step of solving the polynomial equation, how should we deal with the symbol λ while determining the 
//determinant of the matrix.
NUMERICLIBRARY_API void FindEigenValuesAndEigenVectors(double **matrix, int rows, int columns, double *eigenValues, double *eigenVector) {

}


//A=LU
//creates lower diagonal and upper diagonal matrices. The diagonal of the lower diagonal matrice is set to 1.
//order of the algo for 3X3 matrix is:
//iteration 1:  solve for first row of U
//iteration 2:  solve for first column of L
//iteration 3:  solve for second row of U
//iteration 4:  solve for second column of L
//iteration 5:  solve for third row of U
NUMERICLIBRARY_API void LUFactorize(double **matrix, int rows, int columns, double **L, double **U) {
	if (rows != columns) {
		throw std::invalid_argument("Only square matrices can be factorized using LU factorization.");
	}

	//double initCoeffArrL[9] = {
	//	1,	0,	0,
	//	NAN,1,	0,
	//	NAN,NAN,1 };
	double *initCoeffArrL = (double *)malloc((rows*columns) * sizeof(double));

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i == j)
				initCoeffArrL[columns*i + j] = 1;
			else if (i > j)
				initCoeffArrL[columns*i + j] = NAN;
			else if (j > i)
				initCoeffArrL[columns*i + j] = 0;
		}
	}
	Init2DimensionalMatrix(L, rows, columns, initCoeffArrL, rows*columns);

	//double initCoeffArrU[9] = {
	//	NAN,NAN,NAN,
	//	0,  NAN,NAN,
	//	0,  0,	NAN };
	double *initCoeffArrU = (double *)malloc((rows*columns) * sizeof(double));

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i > j)
				initCoeffArrU[columns*i + j] = 0;
			else if (j >= i)
				initCoeffArrU[columns*i + j] = NAN;
		}
	}
	Init2DimensionalMatrix(U, rows, columns, initCoeffArrU, rows*columns);

	int rowIndexerU = 0;
	int colIndexerL = 0;
	for (int i = 0; i < (rows * 2) - 1; i++) {

		double valueToSubtract = 0.0;
		//compute row of U
		if (i % 2 == 0) {
			for (size_t j = 0; j < columns; j++) {
				if (isnan(U[rowIndexerU][j])) {
					for (size_t rowOrColIterator = 0; rowOrColIterator < rows; rowOrColIterator++) {
						if (rowIndexerU != rowOrColIterator) {
							valueToSubtract = valueToSubtract + L[rowIndexerU][rowOrColIterator] * (isnan(U[rowOrColIterator][j]) ? 0.0 : U[rowOrColIterator][j]);
						}
					}
					U[rowIndexerU][j] = (matrix[rowIndexerU][j] - valueToSubtract) / L[rowIndexerU][rowIndexerU];
					valueToSubtract = 0;
				}
			}
			rowIndexerU++;
		}
		else {//compute column of U
			for (size_t j = 0; j < rows; j++) {
				if (isnan(L[j][colIndexerL])) {
					for (size_t rowOrColIterator = 0; rowOrColIterator < columns; rowOrColIterator++) {
						if (colIndexerL != rowOrColIterator) {
							valueToSubtract = valueToSubtract + (isnan(L[j][rowOrColIterator]) ? 0.0 : L[j][rowOrColIterator]) * U[rowOrColIterator][colIndexerL];
						}
					}
					L[j][colIndexerL] = (matrix[j][colIndexerL] - valueToSubtract) / U[colIndexerL][colIndexerL];
					valueToSubtract = 0;
				}
			}
			colIndexerL++;
		}
	}
}

//AX=B. Find X vector
//A=LU
//LUX=B
//suppose UX=T where T is tempVarVector
//use forward substituion in LT=B to find T.
//plug in the value of T in UX=T to find x using backward substitution
//method is also called crout method. How does it differ from cholesky decomposition method?
//To solve AX=B using Cholesky:
//decompose the matrix into: A=LL', 
//then suppose L'X=T where T is tempVarVector
//use forward substituion in LT=B to find T.
//plug in the value of T in L'X=T to find x using backward substitution.
//compared to LU decomposition method, Cholesky decomposition is around 2 times efficient.
NUMERICLIBRARY_API void LUDecompositionLinearEquationsSolver(double **matrix, int rows, int columns, double *variables, double *rhs) {
	double **matL;
	double **matU;
	matL = Create2DimensionalMatrix(rows, columns);
	matU = Create2DimensionalMatrix(rows, columns);

	LUFactorize(matrix, rows, columns, matL, matU);

	//tempVarVector T
	double *T = (double *)malloc((rows) * sizeof(double));

	//find value for T
	ForwardSubstitution(matL, rows, columns, T, rhs);

	//plug in the value of T in UX=T to calculate the varibale verctor
	BackwardSubstitution(matU, rows, columns, variables, T);

	for (int i = 0; i < rows; i++)
		free(matL[i]);
	free(matL);

	for (int i = 0; i < rows; i++)
		free(matU[i]);
	free(matU);
	free(T);
}

//given a matrix, it will reduce the values in all rows under the matrix[rowIndex][rowIndex] value to 0. 
//so the column below matrix[rowIndex][rowIndex] is made 0.
//if row reduction is done for all rows with TriangularForm set to Upper(which is the default), it makes the matrix a upper diagonal matrix
//rhs would be used if we are solving a system of linear equations using Gauss elimination method
NUMERICLIBRARY_API void rowReduction(double **matrix, int rowIndex, int rows, int columns, TriangularForm tf, bool usePivoting, double *rhs) {

	//to reduce round off errors, move a small number out of pivot position
	//other fancy stuff that could be done is to make a fractional number at the pivot position a whole number by 
	//multiplying with a inverse.
	if (usePivoting) {
		if (matrix[rowIndex][rowIndex] == 0 || matrix[rowIndex][rowIndex] < .5)
		{
			//interchange the rows and that interchange should also refelect in rhs vector as well
			throw NotImplementedException("Pivoting is not implemented yet.");
		}
	}

	if (tf == TriangularForm::Upper || tf == TriangularForm::Diagonal) {
		int numberOfRowsToReduce = rows - 1 - rowIndex;

		double rowIndexDataValue = matrix[rowIndex][rowIndex];
		double rowIndexRHSValue;
		if (rhs != nullptr) {
			rowIndexRHSValue = rhs[rowIndex];
		}

		for (size_t i = 1; i <= numberOfRowsToReduce; i++) {
			double rowValueToReduce = matrix[rowIndex + i][rowIndex];
			double rowRHSValueToReduce;
			if (rhs != nullptr) {
				rowRHSValueToReduce = rhs[rowIndex + i];
			}

			double multiplicant = 0;
			if (rowIndexDataValue != 0)
				multiplicant = rowValueToReduce / rowIndexDataValue;


			//multiply the multiplicant rowindex row values and substract it from row to be reduced
			for (size_t j = rowIndex; j < columns; j++) {
				matrix[rowIndex + i][j] = matrix[rowIndex + i][j] - multiplicant * matrix[rowIndex][j];
			}

			//do the same operation on the rhs values as well
			if (rhs != nullptr) {
				rhs[rowIndex + i] = rhs[rowIndex + i] - multiplicant * rhs[rowIndex];
			}
		}
	}

	if (tf == TriangularForm::Lower || tf == TriangularForm::Diagonal) {
		int numberOfRowsToReduce = rowIndex;

		double rowIndexDataValue = matrix[rowIndex][rowIndex];
		double rowIndexRHSValue;
		if (rhs != nullptr) {
			rowIndexRHSValue = rhs[rowIndex];
		}

		for (size_t i = 1; i <= numberOfRowsToReduce; i++) {
			double rowValueToReduce = matrix[rowIndex - i][rowIndex];
			double rowRHSValueToReduce;
			if (rhs != nullptr) {
				rowRHSValueToReduce = rhs[rowIndex - i];
			}

			double multiplicant = 0;
			if (rowIndexDataValue != 0)
				multiplicant = rowValueToReduce / rowIndexDataValue;


			//multiply the multiplicant rowindex row values and substract it from row to be reduced
			for (size_t j = rowIndex; j < columns; j++) {
				matrix[rowIndex - i][j] = matrix[rowIndex - i][j] - multiplicant * matrix[rowIndex][j];
			}

			//do the same operation on the rhs values as well
			if (rhs != nullptr) {
				rhs[rowIndex - i] = rhs[rowIndex - i] - multiplicant * rhs[rowIndex];
			}
		}
	}
	////debug
	//printf("\nRow reduced matrix is:\n");
	//printMatrix(matrix, rows, columns);
	//printf("\n");
}

//rhs would be used if we are solving a system of linear equations using Gauss elimination method which 
//requires the system to be acted upon by row reduction to bring the lhs matrix in upper triangular form
NUMERICLIBRARY_API void reduceToTriangularForm(double **matrix, int rows, int columns, TriangularForm tf, bool usePivoting, double *rhs) {
	for (size_t i = 0; i < rows; i++) {
		rowReduction(matrix, i, rows, columns, tf, usePivoting, rhs);
	}
}

//Gauss elimination method for soving linear equations consists of 2 steps:
//1) read the matrix as well as rhs of the equations
//2)reduce the matrix ot upper triangular form. Pivoting can be used to reduce the round off error by removing 
//  a very small/fractional number from the pivotal position and bringing in a bigger number by either doing a row interchange(partial pivot)
//  or column interchange(complete pivot). Partial pivot requires the same interchange to be done on the rhs and complete pivot requires the 
//  the order of the variables changed in accordance with the columns being interchanged.
//3) use backward substitution to get the solution
NUMERICLIBRARY_API void GaussElimination(double **matrix, int rows, int columns, double *variables, bool usePivoting, double *rhs) {
	if (rows != columns) {
		throw std::invalid_argument("Only square matrices can be evaluated using gauss elimination method.");
	}

	reduceToTriangularForm(matrix, rows, columns, TriangularForm::Upper, usePivoting, rhs);
	BackwardSubstitution(matrix, rows, columns, variables, rhs);
}

//Gauss elimination method for soving linear equations consists of 2 steps:
//1)Read the matrix as well as rhs of the equations
//2)Reduce the matrix to diagonal matrix form(both upper and lower triangular) by using row reduction. You might also do
//  partial and complete pivoting. pivoting is done to reduce round off error. pivot is the 
//  first non-zero element from the left, also called the leading coefficient and pivoting 
//  should be perfromed if the pivot element is a very small number.Pivoting is performed before each row reduction.
//  so if row reduction is eing perfromed for rowindex=0, then before proceeding with row reduction we have to make sure the
//  element a[0][0], the pivot, is not 0 or a very small number so as to avoid round off errors.
//  i think pivoting could also have been performed in gauss elimination method to reduce round off errors.
//  partial pivoting, interchanging rows, affects the rhs as we have to do interchanges in the rhs vector. 
//  Complete pivoting, column interchange, affect the order of variables as we have to intercahnges in the variables vector.
//3)Now here I have reduced the diagonal matrix to identity matrix to get the solution but since the non-diagonal entries are already
//  zero, so using backward substitution or forward substrituion should not amount to much overhead for small matrices to get the solution. 
//  But for large matrices, the nested loops of backward and forward substituion might add a bit of an overhead even if they are just manipulating elements of the matrix that are zero,
//  so reducing the matrix to identity would be a bit faster.
//Note: The total number of calculations required for gauss jordan method are more than gauss elimination. So for large
//      matrices, use gauss elimination
NUMERICLIBRARY_API void GaussJordanElimination(double **matrix, int rows, int columns, double *variables, bool usePivoting, double *rhs) {
	if (rows != columns) {
		throw std::invalid_argument("Only square matrices can be evaluated using gauss elimination method.");
	}

	reduceToTriangularForm(matrix, rows, columns, TriangularForm::Diagonal, usePivoting, rhs);
	reduceToIdentityMatrix(matrix, rows, columns, variables, rhs);

	//they will give the same result with more number of computations. 
	//By reducing the diagonal matrix to identity matrix is computationally cheaper.
	//BackwardSubstitution(matrix, rows, columns, variables, rhs);
	//ForwardSubstitution(matrix, rows, columns, variables, rhs);
}

NUMERICLIBRARY_API void reduceToIdentityMatrix(double **matrix, int rows, int columns, double *variables, double *rhs) {
	if (rows != columns) {
		throw std::invalid_argument("Only square matrices can be reduced to identity matrix.");
	}
	if (!isDiagonalMatrix(matrix, rows, columns)) {
		throw std::invalid_argument("Only diagnonal matrices can be reduced to identity matrix. Use the method reduceToTriangularForm to reduce to diagonal form and then call reduceToIdentityMatrix.");
	}

	for (size_t i = 0; i < rows; i++) {
		if (matrix[i][i] != 1) {
			double inverse = pow(matrix[i][i], -1);
			matrix[i][i] = matrix[i][i] * inverse;
			rhs[i] = rhs[i] * inverse;
			variables[i] = rhs[i];
		}
	}
}

NUMERICLIBRARY_API void rowReduction(double **matrix, int rowIndex, int rows, int columns, TriangularForm tf, bool usePivoting, double **rhs) {
	//to reduce round off errors, move a small number out of pivot position
	//other fancy stuff that could be done is to make a fractional number at the pivot position a whole number by 
	//multiplying with a inverse.
	if (usePivoting) {
		if (matrix[rowIndex][rowIndex] == 0 || matrix[rowIndex][rowIndex] < .5)
		{
			//interchange the rows and that interchange should also refelect in rhs vector as well
			throw NotImplementedException("Pivoting is not implemented yet.");
		}
	}

	if (tf == TriangularForm::Upper || tf == TriangularForm::Diagonal) {
		int numberOfRowsToReduce = rows - 1 - rowIndex;

		double rowIndexDataValue = matrix[rowIndex][rowIndex];
		double rowIndexRHSValue;
		if (rhs != nullptr) {
			rowIndexRHSValue = rhs[rowIndex][rowIndex];
		}

		for (size_t i = 1; i <= numberOfRowsToReduce; i++) {
			double rowValueToReduce = matrix[rowIndex + i][rowIndex];
			double rowRHSValueToReduce;
			if (rhs != nullptr) {
				rowRHSValueToReduce = rhs[rowIndex + i][rowIndex];
			}

			double multiplicant = 0;
			if (rowIndexDataValue != 0)
				multiplicant = rowValueToReduce / rowIndexDataValue;


			//multiply the multiplicant rowindex row values and substract it from row to be reduced
			for (size_t j = rowIndex; j < columns; j++) {
				matrix[rowIndex + i][j] = matrix[rowIndex + i][j] - multiplicant * matrix[rowIndex][j];
			}

			//do the same operation on the rhs values as well
			if (rhs != nullptr) {
				for (size_t j = 0; j < columns; j++) {
					rhs[rowIndex + i][j] = rhs[rowIndex + i][j] - multiplicant * rhs[rowIndex][j];
				}
			}
		}
	}

	if (tf == TriangularForm::Lower || tf == TriangularForm::Diagonal) {
		int numberOfRowsToReduce = rowIndex;

		double rowIndexDataValue = matrix[rowIndex][rowIndex];
		double rowIndexRHSValue;
		if (rhs != nullptr) {
			rowIndexRHSValue = rhs[rowIndex][rowIndex];
		}

		for (size_t i = 1; i <= numberOfRowsToReduce; i++) {
			double rowValueToReduce = matrix[rowIndex - i][rowIndex];
			double rowRHSValueToReduce;
			if (rhs != nullptr) {
				rowRHSValueToReduce = rhs[rowIndex - i][rowIndex];
			}

			double multiplicant = 0;
			if (rowIndexDataValue != 0)
				multiplicant = rowValueToReduce / rowIndexDataValue;


			//multiply the multiplicant rowindex row values and substract it from row to be reduced
			for (size_t j = rowIndex; j < columns; j++) {
				matrix[rowIndex - i][j] = matrix[rowIndex - i][j] - multiplicant * matrix[rowIndex][j];
			}

			//do the same operation on the rhs values as well
			if (rhs != nullptr) {
				for (size_t j = 0; j < columns; j++) {
					rhs[rowIndex - i][j] = rhs[rowIndex - i][j] - multiplicant * rhs[rowIndex][j];
				}
			}
		}
	}
	////debug
	//printf("\nRow reduced matrix is:\n");
	//printMatrix(matrix, rows, columns);
	//printf("\n");
}

//rhs would be used if we are solving a system of linear equations using Gauss elimination method which 
//requires the system to be acted upon by row reduction to bring the lhs matrix in upper triangular form
NUMERICLIBRARY_API void reduceToTriangularForm(double **matrix, int rows, int columns, TriangularForm tf, bool usePivoting, double **rhs) {
	for (size_t i = 0; i < rows; i++) {
		rowReduction(matrix, i, rows, columns, tf, usePivoting, rhs);
	}
}

//1. AA**-1=I
//2. use backward substituion with a slight change.
//3. For rhs u have identity matrix instead of a vector and similarly, instead of the varaibles vector that we have that we
//   have to solve the values for, we have a matrix. The number of columns of the variables matrix give us the number of 
//	 system of simultaneous equations we have to solve for.
NUMERICLIBRARY_API double** InverseOfMatrixUsingGaussElimination(double **matrix, int rows, int columns) {
	if (rows != columns) {
		throw std::invalid_argument("Inverse of only square matrices can be computed using gauss elimination method.");
	}

	double **identityMatrixRHS = Create2DimensionalMatrix(rows, columns);
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i == j)
				identityMatrixRHS[i][j] = 1.0;
			else
				identityMatrixRHS[i][j] = 0.0;
		}
	}

	double **inverseOfMatrix = Create2DimensionalMatrix(rows, columns);
	reduceToTriangularForm(matrix, rows, columns, TriangularForm::Upper, false, identityMatrixRHS);
	BackwardSubstitution(matrix, rows, columns, inverseOfMatrix, identityMatrixRHS);

	for (int i = 0; i < rows; i++)
		free(identityMatrixRHS[i]);
	free(identityMatrixRHS);

	return inverseOfMatrix;
}


//a method of calculating determinants. Much faster than naive determinant evaluation but works only for square matrices
NUMERICLIBRARY_API double PivotalCondensationMethod(double **matrix, int rows, int columns) {
	if (rows != columns) {
		throw std::invalid_argument("Only square matrices can be evaluated using gauss elimination method.");
	}
	if (rows == 1 && columns == 1) {
		return matrix[0][0];
	}
	double v = matrix[0][0];
	rowReduction(matrix, 0, rows, columns, TriangularForm::Upper, false);
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
NUMERICLIBRARY_API double** viewOfMatrix(double **matrix, int rows, int columns, int col) {
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
NUMERICLIBRARY_API double Determinant(double **matrix, int rows, int columns) {
	if (rows != columns) {
		throw std::invalid_argument("Determinant can be calculated only for square matices.");
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

//A diagonal matrix is one in which the main diagonal values are non-zero. Main diagonal is where rowNumber==colNumber
bool isDiagonalMatrix(double **matrix, int rows, int columns) {

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i == j && matrix[i][j] == 0)
				return false;
			else if (i != j && matrix[i][j] != 0)
				return false;
		}
	}
	return true;
}

bool isLowerTriangularMatrix(double **matrix, int rows, int columns) {

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i < j && matrix[i][j] != 0)
				return false;
			//else if (i >= j && matrix[i][j] == 0)
			//	return false;
		}
	}
	return true;
}

bool isUpperTriangularMatrix(double **matrix, int rows, int columns) {

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i > j && matrix[i][j] != 0)
				return false;
			//else if (i <= j && matrix[i][j] == 0)
			//	return false;
		}
	}
	return true;
}

void ForwardSubstitution(double **matrix, int rows, int columns, double *variables, double *rhs) {
	if (!isLowerTriangularMatrix(matrix, rows, columns)) {
		throw std::invalid_argument("Forward substitution only works on lower triangular matrices.");
	}

	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < columns; j++) {
			if (i >= j) {
				double sum = 0;
				for (size_t k = 0; k < i; k++) {
					sum = sum + matrix[i][k] * variables[k];
				}
				variables[i] = (rhs[i] - sum) / matrix[i][i];
				break;
			}
			else
				break;
		}
	}
}

void BackwardSubstitution(double **matrix, int rows, int columns, double *variables, double *rhs) {
	if (!isUpperTriangularMatrix(matrix, rows, columns)) {
		throw std::invalid_argument("Backward substitution works only on upper triangular matrices.");
	}

	for (int i = rows - 1; i >= 0; i--) {
		for (int j = columns - 1; j >= 0; j--) {
			if (i <= j) {
				double sum = 0;
				for (int k = rows - 1; k > i; k--) {
					sum = sum + matrix[i][k] * variables[k];
				}
				variables[i] = (rhs[i] - sum) / matrix[i][i];
				break;
			}
			else
				break;
		}
	}
}

void BackwardSubstitution(double **matrix, int rows, int columns, double **inverseOfMatrix, double **rhsIdentityMatrix) {
	if (!isUpperTriangularMatrix(matrix, rows, columns)) {
		throw std::invalid_argument("Backward substitution works only on upper triangular matrices.");
	}

	for (int l = 0; l < columns; l++) {
		for (int i = rows - 1; i >= 0; i--) {
			for (int j = columns - 1; j >= 0; j--) {
				if (i <= j) {
					double sum = 0;
					for (int k = rows - 1; k > i; k--) {
						sum = sum + matrix[i][k] * inverseOfMatrix[k][l];
					}
					inverseOfMatrix[i][l] = (rhsIdentityMatrix[i][l] - sum) / matrix[i][i];
					break;
				}
				else
					break;
			}
		}
	}
}