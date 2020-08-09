#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\NumericLibrary\NumericLibrary.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

//unittesting as explained here:https://msdn.microsoft.com/en-us/library/hh598953.aspx
//Writing Unit tests for C/C++ with the Microsoft Unit Testing Framework for C++

namespace NumericLibraryTest
{
	TEST_CLASS(InterpolationTests)
	{
	public:

		TEST_METHOD(LinearInterpolationTest)
		{
			//	x		y
			//	---------
			//	1		4
			//	2		5
			//	4		7
			int tableRows = 3;
			int	tableColumns = 2;
			double **interpolationTable = Create2DimensionalMatrix(tableRows, tableColumns);

			double initArr[6] = {
				1,4,
				2,5,
				4,7 };
			Init2DimensionalMatrix(interpolationTable, tableRows, tableColumns, initArr, sizeof(initArr) / sizeof(double));

			double y = LinearInterpolation(3, interpolationTable, 3);

			Assert::AreEqual(6.0, y);

			for (int i = 0; i < tableRows; i++)
				free(interpolationTable[i]);
			free(interpolationTable);
		}

		TEST_METHOD(NewtonCentralDifferencesPolynomialInterpolationTest)
		{
			//odd rows case
			//central difference table for x=5. 
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7											3					0	
			//	1	|4	-3		13									1					0
			//	2	|6	6		28		43							1					1
			//	3	|8	25		56									3					0
			//	4	|10	62											5					0
			//

			wchar_t message[200];

			int intepolationTblRows = 5;
			int	intepolationTblColumns = 2;
			double initArr[10] = {
				2,-7,
				4,-3,
				6,6,
				8,25,
				10,62 };
			double **interpolationTable = Create2DimensionalMatrix(intepolationTblRows, intepolationTblColumns);
			Init2DimensionalMatrix(interpolationTable, intepolationTblRows, intepolationTblColumns, initArr, sizeof(initArr) / sizeof(double));

			double x = 5;
			double y = NewtonPolynomialInterpolation(x, interpolationTable, intepolationTblRows);
			Assert::AreEqual(8.125, y, .001);

			for (int i = 0; i < intepolationTblRows; i++)
				free(interpolationTable[i]);
			free(interpolationTable);

			//even rows case
			//central difference table for x=7.
			//		__________________________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	Δ^5f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7														5					0	
			//	1	|4	-3		13												3					0
			//	2	|6	6		28		43										1					0
			//	3	|8	25		56		27										1					1
			//	4	|10	62		55												3					0
			//	5	|12	80														5					0
			//

			intepolationTblRows = 6;
			intepolationTblColumns = 2;
			double initArr1[12] = {
				2,-7,
				4,-3,
				6,6,
				8,25,
				10,62,
				12,80 };
			interpolationTable = Create2DimensionalMatrix(intepolationTblRows, intepolationTblColumns);
			Init2DimensionalMatrix(interpolationTable, intepolationTblRows, intepolationTblColumns, initArr1, sizeof(initArr1) / sizeof(double));

			x = 7;
			y = NewtonPolynomialInterpolation(x, interpolationTable, intepolationTblRows);
			Assert::AreEqual(7.125, y, .001);

			//even rows case
			//central difference table for x=5.
			//		__________________________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	Δ^5f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7														3					0	
			//	1	|4	-3		13												1					0
			//	2	|6	6		28		43										1					1
			//	3	|8	25		56		27										3					0
			//	4	|10	62		55												5					0
			//	5	|12	80														7					0
			//

			Init2DimensionalMatrix(interpolationTable, intepolationTblRows, intepolationTblColumns, initArr1, sizeof(initArr1) / sizeof(double));

			x = 5;
			y = NewtonPolynomialInterpolation(x, interpolationTable, intepolationTblRows);
			Assert::AreEqual(8.125, y, .001);

			for (int i = 0; i < intepolationTblRows; i++)
				free(interpolationTable[i]);
			free(interpolationTable);
		}

		TEST_METHOD(NewtonPolynomialInterpolationTest)
		{
			//The method describe here is newton difference method
			//If the x points in interpolationTable are equally spaced, then we can use the forward(or backward) difference table
			//for x = 3, the difference table is forward difference as x is near the topmost xi value in the table (3-2 < 10-3 and 3-2 < 6-3).

			//forward difference table for x=3
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7		4		5		5		3			1					1
			//	1	|4	-3		9		10		8					1					0
			//	2	|6	6		19		18							3					0
			//	3	|8	25		37									5					0
			//	4	|10	62											7					0
			//	

			//backward difference table for x=9
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7											7					0
			//	1	|4	-3		4									5					0
			//	2	|6	6		9		5							3					0
			//	3	|8	25		19		10		5					1					0
			//	4	|10	62		37		18		8		3			1					1
			//	

			//central difference table for x=5. 
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7											3					0	
			//	1	|4	-3		13									1					0
			//	2	|6	6		28		43							1					1
			//	3	|8	25		56									3					0
			//	4	|10	62											5					0
			//

			//int tableRows = 5;
			//int	tableColumns = 2;
			//double initArr[10] = {
			//	2,-7,
			//	4,-3,
			//	6,6,
			//	8,25,
			//	10,62 };

			int tableRows = 7;
			int	tableColumns = 2;
			double **interpolationTable = Create2DimensionalMatrix(tableRows, tableColumns);

			double initArr[14] = {
				100,10.63,
				150,13.03,
				200,15.04,
				250,16.81,
				300,18.42,
				350,19.90,
				400,21.27 };
			Init2DimensionalMatrix(interpolationTable, tableRows, tableColumns, initArr, sizeof(initArr) / sizeof(double));

			//this would have used the central difference method but I forced it to use forward differences
			double x = 218;
			double y = NewtonPolynomialInterpolation(x, interpolationTable, tableRows, DifferenceTableDirection::Forward);
			Assert::AreEqual(15.697, y, .001);

			Init2DimensionalMatrix(interpolationTable, tableRows, tableColumns, initArr, sizeof(initArr) / sizeof(double));

			//this used backward differences for polynomial approximation
			x = 410;
			y = NewtonPolynomialInterpolation(x, interpolationTable, tableRows);
			Assert::AreEqual(21.535, y, .001);


			for (int i = 0; i < tableRows; i++)
				free(interpolationTable[i]);
			free(interpolationTable);
		}

		TEST_METHOD(DifferenceTableTest) {
			//forward difference table for x=3
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7		4		5		5		3			1					1
			//	1	|4	-3		9		10		8					1					0
			//	2	|6	6		19		18							3					0
			//	3	|8	25		37									5					0
			//	4	|10	62											7					0
			//	


			wchar_t message[200];

			int intepolationTblRows = 5;
			int	intepolationTblColumns = 2;
			double initArr[10] = {
				2,-7,
				4,-3,
				6,6,
				8,25,
				10,62 };
			double **interpolationTable = Create2DimensionalMatrix(intepolationTblRows, intepolationTblColumns);
			Init2DimensionalMatrix(interpolationTable, intepolationTblRows, intepolationTblColumns, initArr, sizeof(initArr) / sizeof(double));
			int differenceTblCols = 2 + intepolationTblRows + 1;

			//test 1
			double x = 3;
			double **differenceTable = Create2DimensionalMatrix(intepolationTblRows, differenceTblCols);

			DifferenceTableDirection direction = ComputeDifferenceTable(interpolationTable, differenceTable, intepolationTblRows,
				differenceTblCols, x);


			if (DifferenceTableDirection::Forward == direction) {
				Assert::IsTrue(true);
			}
			else {
				_swprintf(message, L"The difference table direction is not forward for %g", x);
				Assert::Fail(message, LINE_INFO());
			}

			//the selected x0 is the first row, x0=2.
			Assert::AreEqual(differenceTable[0][7], 1.0);
			//check the values in the table
			Assert::AreEqual(4.0, differenceTable[0][2]);
			Assert::AreEqual(9.0, differenceTable[1][2]);
			Assert::AreEqual(19.0, differenceTable[2][2]);
			Assert::AreEqual(37.0, differenceTable[3][2]);
			Assert::AreEqual(5.0, differenceTable[0][3]);
			Assert::AreEqual(10.0, differenceTable[1][3]);
			Assert::AreEqual(18.0, differenceTable[2][3]);

			for (int i = 0; i < intepolationTblRows; i++)
				free(differenceTable[i]);
			free(differenceTable);

			//test 2
			//backward difference table for x=9
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7											7					0
			//	1	|4	-3		4									5					0
			//	2	|6	6		9		5							3					0
			//	3	|8	25		19		10		5					1					0
			//	4	|10	62		37		18		8		3			1					1
			//	

			x = 9;
			differenceTable = Create2DimensionalMatrix(intepolationTblRows, differenceTblCols);

			direction = ComputeDifferenceTable(interpolationTable, differenceTable, intepolationTblRows,
				differenceTblCols, x);

			if (DifferenceTableDirection::Backward == direction) {
				Assert::IsTrue(true);
			}
			else {
				_swprintf(message, L"The difference table direction is not backward for %g", x);
				Assert::Fail(message, LINE_INFO());
			}

			//the selected x0 is the last row, x0=4.
			Assert::AreEqual(differenceTable[4][7], 1.0);
			//check the values in the table
			Assert::AreEqual(4.0, differenceTable[1][2]);
			Assert::AreEqual(9.0, differenceTable[2][2]);
			Assert::AreEqual(19.0, differenceTable[3][2]);
			Assert::AreEqual(37.0, differenceTable[4][2]);
			Assert::AreEqual(18.0, differenceTable[4][3]);
			Assert::AreEqual(8.0, differenceTable[4][4]);
			Assert::AreEqual(3.0, differenceTable[4][5]);

			for (int i = 0; i < intepolationTblRows; i++)
				free(differenceTable[i]);
			free(differenceTable);

			//test 3
			//central difference table for x=5. 
			//		______________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7											3					0	
			//	1	|4	-3		13									1					0
			//	2	|6	6		28		43							1					1
			//	3	|8	25		56									3					0
			//	4	|10	62											5					0
			//
			x = 5;
			differenceTable = Create2DimensionalMatrix(intepolationTblRows, differenceTblCols);

			direction = ComputeDifferenceTable(interpolationTable, differenceTable, intepolationTblRows,
				differenceTblCols, x);

			if (DifferenceTableDirection::Center == direction) {
				Assert::IsTrue(true);
			}
			else {
				_swprintf(message, L"The difference table direction is not center for %g", x);
				Assert::Fail(message, LINE_INFO());
			}

			//the selected x0 is the center row, x0=6.
			Assert::AreEqual(differenceTable[2][7], 1.0);
			//check the values in the table
			Assert::AreEqual(differenceTable[1][2], 13.0);
			Assert::AreEqual(differenceTable[2][2], 28.0);
			Assert::AreEqual(differenceTable[3][2], 56.0);
			Assert::AreEqual(differenceTable[2][3], 43.0);

			for (int i = 0; i < intepolationTblRows; i++)
				free(differenceTable[i]);
			free(differenceTable);

			for (int i = 0; i < intepolationTblRows; i++)
				free(interpolationTable[i]);
			free(interpolationTable);
		}

		TEST_METHOD(CentralDifferenceTableTest) {

			wchar_t message[200];

			int intepolationTblRows = 6;
			int intepolationTblColumns = 2;
			double initArr[12] = {
				2,-7,
				4,-3,
				6,6,
				8,25,
				10,62,
				12,80 };
			double **interpolationTable = Create2DimensionalMatrix(intepolationTblRows, intepolationTblColumns);
			Init2DimensionalMatrix(interpolationTable, intepolationTblRows, intepolationTblColumns, initArr, sizeof(initArr) / sizeof(double));
			int differenceTblCols = 2 + intepolationTblRows + 1;

			//test 1. 
			//central difference table for x=7.
			//		__________________________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	Δ^5f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7														5					0	
			//	1	|4	-3		13												3					0
			//	2	|6	6		28		43										1					0
			//	3	|8	25		56		27										1					1
			//	4	|10	62		55												3					0
			//	5	|12	80														5					0
			//
			int x = 7;
			double **differenceTable = Create2DimensionalMatrix(intepolationTblRows, differenceTblCols);

			DifferenceTableDirection direction = ComputeDifferenceTable(interpolationTable, differenceTable, intepolationTblRows,
				differenceTblCols, x);

			if (DifferenceTableDirection::Center == direction) {
				Assert::IsTrue(true);
			}
			else {
				_swprintf(message, L"The difference table direction is not center for %g", x);
				Assert::Fail(message, LINE_INFO());
			}

			//the selected x0 is the center row, x0=6.
			Assert::AreEqual(differenceTable[3][8], 1.0);
			//check values in the table
			Assert::AreEqual(13.0, differenceTable[1][2]);
			Assert::AreEqual(28.0, differenceTable[2][2]);
			Assert::AreEqual(56.0, differenceTable[3][2]);
			Assert::AreEqual(55.0, differenceTable[4][2]);
			Assert::AreEqual(43.0, differenceTable[2][3]);
			Assert::AreEqual(27.0, differenceTable[3][3]);

			for (int i = 0; i < intepolationTblRows; i++)
				free(differenceTable[i]);
			free(differenceTable);

			//test 2. 
			//central difference table for x=5.
			//		__________________________________________________________________________________________
			//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	Δ^5f(xi)	DistanceOfxiFromX	Selectedx0
			//	0	|2	-7														3					0	
			//	1	|4	-3		13												1					0
			//	2	|6	6		28		43										1					1
			//	3	|8	25		56		27										3					0
			//	4	|10	62		55												5					0
			//	5	|12	80														7					0
			//

			x = 5;
			differenceTable = Create2DimensionalMatrix(intepolationTblRows, differenceTblCols);

			direction = ComputeDifferenceTable(interpolationTable, differenceTable, intepolationTblRows,
				differenceTblCols, x);

			if (DifferenceTableDirection::Center == direction) {
				Assert::IsTrue(true);
			}
			else {
				_swprintf(message, L"The difference table direction is not center for %g", x);
				Assert::Fail(message, LINE_INFO());
			}

			//the selected x0 is the center row, x0=6.
			Assert::AreEqual(differenceTable[2][8], 1.0);
			//check values in the table
			Assert::AreEqual(13.0, differenceTable[1][2]);
			Assert::AreEqual(28.0, differenceTable[2][2]);
			Assert::AreEqual(56.0, differenceTable[3][2]);
			Assert::AreEqual(55.0, differenceTable[4][2]);
			Assert::AreEqual(43.0, differenceTable[2][3]);
			Assert::AreEqual(27.0, differenceTable[3][3]);

			for (int i = 0; i < intepolationTblRows; i++)
				free(differenceTable[i]);
			free(differenceTable);


			for (int i = 0; i < intepolationTblRows; i++)
				free(interpolationTable[i]);
			free(interpolationTable);
		}
	};
}