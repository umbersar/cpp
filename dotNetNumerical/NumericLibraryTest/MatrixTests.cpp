#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\NumericLibrary\NumericLibrary.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

//unittesting as explained here:https://msdn.microsoft.com/en-us/library/hh598953.aspx
//Writing Unit tests for C/C++ with the Microsoft Unit Testing Framework for C++

namespace NumericLibraryTest
{
	TEST_CLASS(MatrixTests)
	{
	public:

		TEST_METHOD(IsDiagonalMatrixTest)
		{
			int matrixRows = 4;
			int	matrixColumns = 4;
			double **matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initArr[16] = {
				2,0,0,0,
				0,1,0,0,
				0,0,1,0,
				0,0,0,1 };
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, initArr, sizeof(initArr) / sizeof(double));

			Assert::AreEqual(true, isDiagonalMatrix(matA, matrixRows, matrixColumns));

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
		}

		TEST_METHOD(isLowerTriangularMatrixTest)
		{
			int matrixRows = 4;
			int	matrixColumns = 4;
			double **matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initArr[16] = {
				2,0,0,0,
				2,1,0,0,
				2,2,1,0,
				2,2,2,1 };
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, initArr, sizeof(initArr) / sizeof(double));

			Assert::AreEqual(true, isLowerTriangularMatrix(matA, matrixRows, matrixColumns));

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
		}

		TEST_METHOD(isUpperTriangularMatrixTest)
		{
			int matrixRows = 4;
			int	matrixColumns = 4;
			double **matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initArr[16] = {
				2,1,2,3,
				0,1,2,3,
				0,0,1,2,
				0,0,0,1 };
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, initArr, sizeof(initArr) / sizeof(double));

			Assert::AreEqual(true, isUpperTriangularMatrix(matA, matrixRows, matrixColumns));

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
		}

		TEST_METHOD(ForwardSubstitutionLinearAlgebraTest)
		{
			int matrixRows = 4;
			int	matrixColumns = 4;
			double **matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[16] = {
				3,0,0,0,
				-1,1,0,0,
				3,-2,-1,0,
				1,-2,6,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables[4];
			double rhs[4] = { 5,6,4,2 };
			ForwardSubstitution(matCoefficients, matrixRows, matrixColumns, variables, rhs);

			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);

			//todo: can u do asserts on full array
			Assert::AreEqual(1.6666666666666667, variables[0]);
			Assert::AreEqual(7.6666666666666670, variables[1]);
			Assert::AreEqual(-14.333333333333334, variables[2]);
			Assert::AreEqual(50.833333333333336, variables[3]);
		}

		TEST_METHOD(BackwardSubstitutionLinearAlgebraTest)
		{
			int matrixRows = 3;
			int	matrixColumns = 3;
			double **matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				1,-2,1,
				0,1,6,
				0,0,1 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables[3];
			double rhs[3] = { 4,-1,2 };
			BackwardSubstitution(matCoefficients, matrixRows, matrixColumns, variables, rhs);

			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);

			//todo: can u do asserts on full array
			Assert::AreEqual(-24.0, variables[0]);
			Assert::AreEqual(-13.0, variables[1]);
			Assert::AreEqual(2.0, variables[2]);
		}

		TEST_METHOD(DeterminantOfMatrixUsingNaiveMethodTest)
		{
			int matrixRows = 3;
			int matrixColumns = 3;
			double seed = 0;
			double **matA;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

			printf("The input matrix for evaluating determinant is:\n");
			printMatrix(matA, matrixRows, matrixColumns);

			double det = Determinant(matA, matrixRows, matrixColumns);
			printf("The determinant of the matrix using naive determinant method is: %lf\n", det);

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);

			Assert::AreEqual(0.0, det);

			matrixRows = 2;
			matrixColumns = 2;
			seed = 0;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

			printf("The input matrix for evaluating determinant is:\n");
			printMatrix(matA, matrixRows, matrixColumns);

			det = Determinant(matA, matrixRows, matrixColumns);
			printf("The determinant of the matrix using naive determinant method is: %lf\n", det);

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
			printf("\n");

			Assert::AreEqual(-2.0, det);
		}

		TEST_METHOD(RunTimeComparisonsDeterminantOfMatrixUsingNaiveMethodAndEigenLibraryTest)
		{
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matA;

			cpu_timer timer;
			timer.stop();

			//do some runtime comparions between my code and eigen code for big matrix determinant calculations
			matrixRows = 10;
			matrixColumns = 10;
			seed = 0;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

			/*printf("The input matrix for evaluating determinant is:\n");
			printMatrix(matA, matrixRows, matrixColumns);*/
			timer.start();
			printf("The determinant of the matrix using naive determinant method is: %lf\n", Determinant(matA, matrixRows, matrixColumns));
			timer.stop();

			cpu_times timesNaive = timer.elapsed();
			printf("Timings for evaluating determinant of the matrix using my code %s\n", timer.format().c_str());
			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
			printf("\n");

			//eigen library
			//Matrix<double, 10, 10>  matE;
			MatrixXd matE(matrixRows, matrixColumns);

			for (size_t i = 0; i < matrixRows; i++) {
				for (size_t j = 0; j < matrixColumns; j++) {
					matE(i, j) = matrixColumns*i + j;
				}
			}

			timer.start();
			printf("The determinant using Eigen library is: %lf\n", matE.determinant());
			timer.stop();

			cpu_times timesEigen = timer.elapsed();
			printf("Timings for evaluating determinant of the matrix using eigen library is %s\n", timer.format().c_str());
			printf("\n");

			Assert::IsTrue(timesNaive.wall > timesEigen.wall);
		}

		TEST_METHOD(DeterminantOfMatrixUsingPivotalCondensationMethodTest)
		{
			//Determinant of the matrix using pivotal condensation method.
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matA;

			matrixRows = 3;
			matrixColumns = 3;
			seed = 1;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

			printf("The input matrix for evaluating determinant using pivotal condensation method is:\n");
			printMatrix(matA, matrixRows, matrixColumns);
			double det = PivotalCondensationMethod(matA, matrixRows, matrixColumns);
			printf("And the determinant of the matrix using pivotal condensation method is: %lf\n", det);
			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
			printf("\n");

			Assert::AreEqual(0.0, det, .001);

			//another example of pivotal condenstion determinant evaluation
			matrixRows = 4;
			matrixColumns = 4;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initArr[16] = { 2,4,6,8,3,1,2,1,1,2,-1,2,2,3,4,1 };
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, initArr, sizeof(initArr) / sizeof(double));

			printf("The input matrix for evaluating determinant using pivotal condensation method is:\n");
			printMatrix(matA, matrixRows, matrixColumns);
			det = PivotalCondensationMethod(matA, matrixRows, matrixColumns);
			printf("And the determinant of the matrix using pivotal condensation method is: %lf\n", det);
			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
			printf("\n");

			Assert::AreEqual(-180, det, .001);
		}

		TEST_METHOD(RunTimeComaprisonsDeterminantOfMatrixUsingPivotalCondensationMethodAndNaiveMethodTest)
		{
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matA;

			cpu_timer timer;
			timer.stop();

			matrixRows = 10;
			matrixColumns = 10;
			seed = 1;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

			//printf("The input matrix for evaluating determinant using pivotal condensation method is:\n");
			//printMatrix(matA, matrixRows, matrixColumns);

			timer.start();
			double det = PivotalCondensationMethod(matA, matrixRows, matrixColumns);
			timer.stop();
			cpu_times timesPivot = timer.elapsed();
			printf("Timings for evaluating determinant of the matrix using pivotal condensation method is %s", timer.format().c_str());
			printf("And the determinant of the matrix using pivotal condensation method is: %lf\n", det);
			printf("\n");

			//reset the matrix as pivotal condensation method reduces/manipulates the matrix
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);


			/*printf("The input matrix for evaluating determinant using naive determinant method is:\n");
			printMatrix(matA, matrixRows, matrixColumns);*/
			timer.start();
			det = Determinant(matA, matrixRows, matrixColumns);
			timer.stop();
			cpu_times timesNaive = timer.elapsed();
			printf("Timings for evaluating determinant of the matrix using my code %s", timer.format().c_str());
			printf("And the determinant of the matrix using naive determinant method is: %lf\n", det);

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);
			printf("\n");

			Assert::IsTrue(timesNaive.wall > timesPivot.wall);
		}

		TEST_METHOD(GaussEliminationLinearAlgebraTest)
		{
			int matrixRows = 3;
			int	matrixColumns = 3;
			double **matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				2,1,-1,
				-3,-1,2,
				-2,1,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables[3];
			double rhs[3] = { 8,-11,-3 };


			GaussElimination(matCoefficients, matrixRows, matrixColumns, variables, false, rhs);
			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
			printf("\n");

			//todo: can u do asserts on full array
			Assert::AreEqual(2.0, variables[0]);
			Assert::AreEqual(3.0, variables[1]);
			Assert::AreEqual(-1.0, variables[2]);

			//another test case
			matrixRows = 3;
			matrixColumns = 3;
			matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr1[9] = {
				2,-1,2,
				1,-2,1,
				3,-1,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr1, sizeof(initCoeffArr1) / sizeof(double));

			double variables1[3];
			double rhs1[3] = { 10,8,11 };


			GaussElimination(matCoefficients, matrixRows, matrixColumns, variables1, false, rhs1);
			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
			printf("\n");

			//todo: can u do asserts on full array
			Assert::AreEqual(1.0, variables1[0]);
			Assert::AreEqual(-2.0, variables1[1]);
			Assert::AreEqual(3.0, variables1[2]);
		}

		TEST_METHOD(GaussJordanEliminationLinearAlgebraTest)
		{
			int matrixRows = 3;
			int	matrixColumns = 3;
			double **matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				2,1,-1,
				-3,-1,2,
				-2,1,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables[3];
			double rhs[3] = { 8,-11,-3 };


			GaussJordanElimination(matCoefficients, matrixRows, matrixColumns, variables, false, rhs);
			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
			printf("\n");

			//todo: can u do asserts on full array
			Assert::AreEqual(2.0, variables[0]);
			Assert::AreEqual(3.0, variables[1]);
			Assert::AreEqual(-1.0, variables[2]);

			//another test case
			matrixRows = 3;
			matrixColumns = 3;
			matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr1[9] = {
				2,-1,2,
				1,-2,1,
				3,-1,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr1, sizeof(initCoeffArr1) / sizeof(double));

			double variables1[3];
			double rhs1[3] = { 10,8,11 };


			GaussJordanElimination(matCoefficients, matrixRows, matrixColumns, variables1, false, rhs1);
			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
			printf("\n");

			//todo: can u do asserts on full array
			Assert::AreEqual(1.0, variables1[0]);
			Assert::AreEqual(-2.0, variables1[1]);
			Assert::AreEqual(3.0, variables1[2]);
		}

		TEST_METHOD(ExceptionsTest)
		{
			//char buffer[100];
			wchar_t message[200];

			int matrixRows = 3;
			int	matrixColumns = 3;
			double **matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				2,1,-1,
				-3,-1,2,
				-2,1,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables[3];
			double rhs[3] = { 8,-11,-3 };

			try
			{
				GaussJordanElimination(matCoefficients, matrixRows, matrixColumns, variables, true, rhs);
				//sprintf(buffer, "Expected NotImplementedException");
				_swprintf(message, L"Expected NotImplementedException");

				//// Convert to a wchar_t*
				//size_t origsize = strlen(buffer) + 1;
				//const size_t newsize = 100;
				//size_t convertedChars = 0;
				//wchar_t wcstring[newsize];
				//mbstowcs_s(&convertedChars, wcstring, origsize, buffer, _TRUNCATE);
				//wcscat_s(wcstring, L" (wchar_t *)");

				Assert::Fail(message, LINE_INFO());
			}
			catch (NotImplementedException ex)
			{
				//swallow the exception
				//ex.what();
			}
			catch (...)
			{
				//sprintf(buffer, "Incorrect exception.");

				//// Convert to a wchar_t*
				//size_t origsize = strlen(buffer) + 1;
				//const size_t newsize = 100;
				//size_t convertedChars = 0;
				//wchar_t wcstring[newsize];
				//mbstowcs_s(&convertedChars, wcstring, origsize, buffer, _TRUNCATE);
				//wcscat_s(wcstring, L" (wchar_t *)");

				Assert::Fail(message, LINE_INFO());
			}
			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
		}

		TEST_METHOD(CompareGaussJordanAndGaussEliminationRunTimesLinearAlgebraTest)
		{
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matCoefficients;

			cpu_timer timer;
			timer.stop();

			matrixRows = 3;
			matrixColumns = 3;
			matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				2,1,-1,
				-3,-1,2,
				-2,1,2 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables[3];
			double rhs[3] = { 8,-11,-3 };

			timer.start();
			GaussJordanElimination(matCoefficients, matrixRows, matrixColumns, variables, false, rhs);
			timer.stop();
			cpu_times gaussJordanTime = timer.elapsed();

			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
			printf("\n");

			//todo: can u do asserts on full array
			Assert::AreEqual(2.0, variables[0]);
			Assert::AreEqual(3.0, variables[1]);
			Assert::AreEqual(-1.0, variables[2]);

			//another test case
			matrixRows = 3;
			matrixColumns = 3;
			matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			/*double initCoeffArr[9] = {
				2,1,-1,
				-3,-1,2,
				-2,1,2 };*/
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double variables1[3];
			double rhs1[3] = { 8,-11,-3 };

			timer.start();
			GaussElimination(matCoefficients, matrixRows, matrixColumns, variables1, false, rhs1);
			timer.stop();
			cpu_times gaussTime = timer.elapsed();

			for (int i = 0; i < matrixRows; i++)
				free(matCoefficients[i]);
			free(matCoefficients);
			printf("\n");

			//todo: can u do asserts on full array
			Assert::AreEqual(2.0, variables[0]);
			Assert::AreEqual(3.0, variables[1]);
			Assert::AreEqual(-1.0, variables[2]);

			Assert::IsTrue(gaussJordanTime.wall > gaussTime.wall);
		}

		TEST_METHOD(InverseOfMatrixUsingGaussEliminationTest) {
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matCoefficients;

			matrixRows = 3;
			matrixColumns = 3;
			matCoefficients = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				1,2,3,
				2,5,3,
				1,0,8 };
			Init2DimensionalMatrix(matCoefficients, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			double **inverse = InverseOfMatrixUsingGaussElimination(matCoefficients, matrixRows, matrixColumns);

			Assert::AreEqual(-40.0, inverse[0][0]);
			Assert::AreEqual(13.0, inverse[1][0]);
			Assert::AreEqual(5.0, inverse[2][0]);

			Assert::AreEqual(16.0, inverse[0][1]);
			Assert::AreEqual(-5.0, inverse[1][1]);
			Assert::AreEqual(-2.0, inverse[2][1]);

			Assert::AreEqual(9.0, inverse[0][2]);
			Assert::AreEqual(-3.0, inverse[1][2]);
			Assert::AreEqual(-1.0, inverse[2][2]);

			for (int i = 0; i < matrixRows; i++)
				free(inverse[i]);
			free(inverse);
		}

		TEST_METHOD(LUFactorizeTest) {
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matA;
			double **matL;
			double **matU;

			matrixRows = 3;
			matrixColumns = 3;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			matL = Create2DimensionalMatrix(matrixRows, matrixColumns);
			matU = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double initCoeffArr[9] = {
				3,2,7,
				2,3,1,
				3,4,1 };
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			LUFactorize(matA, matrixRows, matrixColumns, matL, matU);

			Assert::AreEqual(1.0, matL[0][0]);
			Assert::AreEqual(2.0 / 3, matL[1][0]);
			Assert::AreEqual(1.0, matL[2][0]);

			Assert::AreEqual(0.0, matL[0][1]);
			Assert::AreEqual(1.0, matL[1][1]);
			Assert::AreEqual(6.0 / 5, matL[2][1]);

			Assert::AreEqual(0.0, matL[0][2]);
			Assert::AreEqual(0.0, matL[1][2]);
			Assert::AreEqual(1.0, matL[2][2]);

			Assert::AreEqual(3.0, matU[0][0]);
			Assert::AreEqual(0.0 / 3, matU[1][0]);
			Assert::AreEqual(0.0, matU[2][0]);

			Assert::AreEqual(2.0, matU[0][1]);
			Assert::AreEqual(5.0 / 3, matU[1][1]);
			Assert::AreEqual(0.0 / 5, matU[2][1]);

			Assert::AreEqual(7.0, matU[0][2]);
			Assert::AreEqual(-11.0 / 3, matU[1][2], .001);
			Assert::AreEqual(-8.0 / 5, matU[2][2], .001);

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);

			for (int i = 0; i < matrixRows; i++)
				free(matL[i]);
			free(matL);

			for (int i = 0; i < matrixRows; i++)
				free(matU[i]);
			free(matU);
		}

		TEST_METHOD(LULinearEquationSolverTest) {
			int matrixRows;
			int matrixColumns;
			double seed;
			double **matA;

			matrixRows = 3;
			matrixColumns = 3;
			matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

			double variables[3];
			double rhs[3] = { 4,5,7 };

			double initCoeffArr[9] = {
				3,2,7,
				2,3,1,
				3,4,1 };
			Init2DimensionalMatrix(matA, matrixRows, matrixColumns, initCoeffArr, sizeof(initCoeffArr) / sizeof(double));

			LUDecompositionLinearEquationsSolver(matA, matrixRows, matrixColumns, variables, rhs);

			Assert::AreEqual(7.0 / 8, variables[0], .00000001);
			Assert::AreEqual(9.0 / 8, variables[1], .00000001);
			Assert::AreEqual(-1.0 / 8, variables[2], .00000001);

			for (int i = 0; i < matrixRows; i++)
				free(matA[i]);
			free(matA);

		}
	};
}