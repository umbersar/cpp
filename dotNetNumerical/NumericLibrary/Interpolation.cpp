#include "stdafx.h"
#include "NumericLibrary.h"

//interpolation: when we have a value table of values of x and y, where y is f(x) but we do not know the f(x), then 
//to calculate a value of y for any other given x is called interpolation or extrapolation.
//	x		y
//	---------
//	1		4
//	2		5
//	4		7
//here, if we have to calculate the value of y for x=3 it is called interpolation and y for x=5 is called extrapolation.
//the word interpolate means estimate. 
//Estimation can be done by regression or curve fitting. In curve fitting problems, the constraint that the
//interpolant(estimate) has to go exactly through the data points is relaxed. It is only required to approach the 
//data points as closely as possible (within some other constraints). This requires parameterizing the potential 
//interpolants and having some way of measuring the error. In the simplest case this leads to least squares approximation.


//taylor series is f(x) = f(x0) + (f'(x0)/1!)(x- x0) + (f''(x0)/2!) (x- x0)^2 + (f'''(0)/3!) (x- x0)^3
//the first 2 terms: f(x) = f(x0) + (f'(x0)/1!)(x- x0). Now f'(x) is slope and since it is linear interpolation, the slope
//will remain same at any point on the curve. So solpe can be calculate using (y1-y0)/(x1-x0) and rewritting the formula it 
//becomes:y = f(x) = y0 + (y1-y0)/(x1-x0) *  (x- x0)
double LinearInterpolation(double x, double **interpolationTable, int intepolationTblRows) {

	if (intepolationTblRows < 2)
		throw std::invalid_argument("Linear interpolation cannot proceed with less than 2 entries in the interpolaion "
			"table as 2 entries are need to calculate the linear slope.");

	//verify that the slope is linear based on the existing entries in interpolationTable
	double slope = NAN;
	for (size_t i = 0; i < intepolationTblRows - 1; i++) {
		//todo:cehck for divide by zero
		double tempSlope = (interpolationTable[i + 1][1] - interpolationTable[i][1]) / (interpolationTable[i + 1][0] - interpolationTable[i][0]);
		if (!isnan(slope)) {
			if (slope != tempSlope) {
				throw std::invalid_argument("InterpolationTable entries do not compute to linear slopes.");
			}
		}
		slope = tempSlope;
	}

	double y = interpolationTable[0][1] + slope *  (x - interpolationTable[0][0]);
	return y;
}

//taylor series is f(x) = f(x0) + (f'(x0)/1!)(x- x0) + (f''(x0)/2!) (x- x0)^2 + (f'''(0)/3!) (x- x0)^3 + ...
//Now f'(x) is slope and solpe can be calculate using (y1-y0)/(x1-x0) and rewritting the formula it 
//becomes:y = f(x) = y0 + (y1-y0)/(x1-x0) *  (x- x0) for first 2 terms. But if we do not use the first 2 terms but 
//instead use the polynomial to do polynomial iterpolation in place of linear inpterpolation, then how do we find f''(x0) 
//and f'''(x0) and so on. Since we do not not the function, we cannot find the derivatives. only first order derivative 
//can be simulated by finding slope. So what to do with higher order dervatives?
//The idea is to replace the derivatives with finite differences that approximate them. The resulting methods are called finite difference methods.
//Now first order difference Δf(xi) for i=0 could be written as y1-y0 and if we divide it by x1-x0, should't we get the first
//order derivative f'(x0). if yes, then what is the relation between higher order derivatives and higher order differences? Maybe
//we will talk about the relationship finite differences and derivatives later. Fow now lets assume given a interpolation table, the
//the nth degree polynomial approximation for n data points is(explained here https://www3.nd.edu/~coast/jjwteach/www/www/30125/pdfnotes/lecture4_7v14.pdf I have also saved it in git NewtonPolynomialInterpolation.pdf):
//f(x) = f(x0) + (Δf(x0)/((x1-x0)*1!))(x- x0) + (Δ²f(x0)/((x1-x0)²*2!)) (x- x0)*(x-x1) + (Δ³f(0)/((x1-x0)³*3!)) (x- x0)*(x-x1)*(x-x2) + ...
//Now we can use the equation for polynomial approximation given above if we do not know the function(hence cannot use taylor series because we have to calculate the derivatives of the f)
//but we are given some x and f(x) results and using them we have to calculate f(x) for a point not in the interpolation table

//Given below is the mapping between difference table and derivatives:
//If the x points in interpolationTable are equally spaced, then we can use the forward(or backward) difference table
//for x = 3, the difference table is forward difference as x is near the topmost xi value in the table (3-2 < 10-3 and 3-2 < 6-3).
//
//forward difference table for x=3
//		______________________________________________________________________________
//	i	|xi	f(xi)	Δf(xi)	Δ²f(xi)	Δ³f(xi)	Δ^4f(xi)	DistanceOfxiFromX	Selectedx0
//	0	|2	-7		4		5		5		3			1					1
//	1	|4	-3		9		10		8					1					0
//	2	|6	6		19		18							3					0
//	3	|8	25		37									5					0
//	4	|10	62											7					0
//f'(x0)=(f(x0 + delta) - f(x0)) / ((x0 + delta) - x0) = Δf(x0) / delta
//f''(x0) = Δ²f(x0)/delta² and so on and can be generalized as:
//fⁿ(x0) = Δⁿf(x0)/deltaⁿ where x0 could be any selected row in the table
//
//backward difference table for x=9
//		______________________________________________________________________________
//	i	|xi	f(xi)	∇f(xi)	∇²f(xi)	∇³f(xi)	∇^4f(xi)	DistanceOfxiFromX	Selectedx0
//	0	|2	-7											7					0
//	1	|4	-3		4									5					0
//	2	|6	6		9		5							3					0
//	3	|8	25		19		10		5					1					0
//	4	|10	62		37		18		8		3			1					1
//	f'(x0)=(f(x0) - f(x0 - delta)) / (x0 - (x0 - delta)) =  ∇f(x0) / delta
//f''(x0) = ∇²f(x0)/delta² and so on and can be generalized as:
//fⁿ(x0) = ∇ⁿf(x0)/deltaⁿ where x0 could be any selected row in the table
//
//central difference table for x=5.
//		______________________________________________________________________________
//	i	|xi	f(xi)	δf(xi)	δ²f(xi)	δ³f(xi)	δ^4f(xi)	DistanceOfxiFromX	Selectedx0
//	0	|2	-7											3					0	
//	1	|4	-3		13									1					0
//	2	|6	6		28		43							1					1
//	3	|8	25		56									3					0
//	4	|10	62											5					0
//	f'(x0)=(f(x0 + delta) - f(x0 - delta)) / ((x0 + delta) - (x0 - delta)) = δf(x0) / 2delta
//f''(x0) = δ²f(x0)/2delta² and so on and can be generalized as:z
//fⁿ(x0) = δⁿf(x0)/2deltaⁿ where x0 could be any selected row in the table
//note that for central difference, f'(x0) could also be expressed in as an airthematic average of forward and backward differences:
//f'(x0)= (1/2) * [(Δf(x0) / delta) + (∇f(x0) / delta)] = (1/2) * ((Δf(x0) + ∇f(x0) ) / delta) = ( Δf(x0) + ∇f(x0) ) / 2*delta
//which could be generalized as for odd or even number of derivative order: 
//1. fⁿ(xi)= ( Δⁿf(x[i+n/2]) + ∇ⁿf(x[i-n/2]) ) / 2*deltaⁿ if the n(derivative order) is even
//2. fⁿ(xi)= ( Δⁿf(x[i+(n-1)/2]) + ∇ⁿf(x[i-(n-1)/2]) ) / 2*deltaⁿ if the n(derivative order) is odd
//So derivatives for central differences can be found using both the central difference table as well as an airthemative average of forward 
//and backward difference
//much of above info gathered from: https://www3.nd.edu/~coast/jjwteach/www/www/30125/pdfnotes/lecture7_12v09.pdf or in git DerivativesFromFiniteDifferences.pdf

double NewtonPolynomialInterpolation(double x, double **interpolationTable, int intepolationTblRows, DifferenceTableDirection diffTblDirection) {

	if (intepolationTblRows < 2) {
		throw std::invalid_argument("Linear interpolation cannot proceed with less than 2 entries in the interpolaion "
			"table as 2 entries are need to calculate the linear slope.");
	}
	//verify that the slope is linear based on the existing entries in interpolationTable
	double slope = NAN;
	bool isSlopeSame = true;
	for (size_t i = 0; i < intepolationTblRows - 1; i++) {
		//todo:cehck for divide by zero
		double tempSlope = (interpolationTable[i + 1][1] - interpolationTable[i][1]) / (interpolationTable[i + 1][0] - interpolationTable[i][0]);
		if (!isnan(slope)) {
			if (slope != tempSlope) {
				isSlopeSame = false;
				break;
			}
		}
		slope = tempSlope;
	}

	if (isSlopeSame == true) {
		throw std::invalid_argument("InterpolationTable entries compute to linear slopes. Hence Linear interpolation would be more efficient");
	}

	double lastSpacing;
	double newSpacing;
	for (size_t i = 0; i < intepolationTblRows - 1; i++) {
		newSpacing = interpolationTable[i + 1][0] - interpolationTable[i][0];
		if (i > 0) {
			if (newSpacing != lastSpacing) {
				throw std::invalid_argument("Spacing for x values in the InterpolationTable is not equal. Polynomial interpolation using newton forward or backward method works only with equally spaced x entries.");
			}
		}
		lastSpacing = newSpacing;
	}

	int differenceTblCols = 2 + intepolationTblRows + 1;// -1;
	double **differenceTable = Create2DimensionalMatrix(intepolationTblRows, differenceTblCols);// -1);

	DifferenceTableDirection direction = ComputeDifferenceTable(interpolationTable, differenceTable, intepolationTblRows, differenceTblCols, x, diffTblDirection);

	double sumOfTerms;
	if (direction == DifferenceTableDirection::Forward) {
		//now implement f(x) = f(x0) + (Δf(x0)/((x1-x0)*1!))(x- x0) + (Δ²f(x0)/((x1-x0)²*2!)) (x- x0)*(x-x1) + (Δ³f(0)/((x1-x0)³*3!)) (x- x0)*(x-x1)*(x-x2) + ...
		//now instead of taking x0 as the last x from the inerpolation table(this is backward differences), use the x value from the inerpolation table nearest to the 'x' we are finding value for
		//So x0 need not be the bottom-most x value in the table. ComputeDifferenceTable has already given you the direction of the difference table, 
		//the row index of x0 value in the table and the backward differences in the row.
		//based on the position of x0 in the table, we also have to adjust the while loop counter below

		//so find the x0 nearest to x from the last column of the DifferenceTable.
		int rowIndexNearestXI = 0;

		for (int i = 0; i < intepolationTblRows; i++) {
			if (differenceTable[i][differenceTblCols - 1] == 1) {
				rowIndexNearestXI = i;
				break;
			}
		}

		int index = 1;
		double newValue;
		double x0 = differenceTable[rowIndexNearestXI][0];
		double x1 = differenceTable[rowIndexNearestXI + 1][0];
		double delta = x1 - x0;
		sumOfTerms = differenceTable[rowIndexNearestXI][1];
		//for the selected row x0, how far right in the row you can go inside the differencetable to use higher order differences.
		while (index < intepolationTblRows - rowIndexNearestXI) {

			double diff = differenceTable[rowIndexNearestXI][1 + index];
			double deltaPower = pow(delta, index);
			double fact = factorial(index);
			double multiplicant;
			double multiplicantContainer;
			for (int i = 0; i < index; i++) {
				multiplicant = x - differenceTable[rowIndexNearestXI + i][0];
				if (i != 0) {
					multiplicant = multiplicant*multiplicantContainer;
				}
				multiplicantContainer = multiplicant;
			}

			newValue = (diff / (deltaPower*fact))*multiplicantContainer;

			sumOfTerms = sumOfTerms + newValue;
			index++;
		}
	}
	else if (direction == DifferenceTableDirection::Backward) {
		//now implement f(x) = f(x0) + (Δf(x0)/((x1-x0)*1!))(x- x0) + (Δ²f(x0)/((x1-x0)²*2!)) (x- x0)*(x-x1) + (Δ³f(0)/((x1-x0)³*3!)) (x- x0)*(x-x1)*(x-x2) + ...
		//now instead of taking x0 as the last x from the inerpolation table(this is backward differences), use the x value from the inerpolation table nearest to the 'x' we are finding value for
		//So x0 need not be the bottom-most x value in the table. ComputeDifferenceTable has already given you the direction of the difference table, 
		//the row index of x0 value in the table and the backward differences in the row.
		//based on the position of x0 in the table, we also have to adjust the while loop counter below

		//texts refer to x0 for backward difference as xn. It does not make a difference as x0 is supposed to be any starting value. 
		//does not make a difference to the formula. The only thing to remember is that we are working backwards in the table.

		//so find the x0 nearest to x from the last column of the DifferenceTable.
		int rowIndexNearestXI = intepolationTblRows - 1;

		for (int i = intepolationTblRows - 1; i >= 0; i++) {
			if (differenceTable[i][differenceTblCols - 1] == 1) {
				rowIndexNearestXI = i;
				break;
			}
		}
		int index = 1;
		double newValue;
		double x0 = differenceTable[rowIndexNearestXI][0];
		double x1 = differenceTable[rowIndexNearestXI - 1][0];
		double delta = x0 - x1;
		sumOfTerms = differenceTable[rowIndexNearestXI][1];
		//for the selected row x0, how far right in the row you can go inside the differencetable to use higher order differences.
		while (index < rowIndexNearestXI + 1) {

			double diff = differenceTable[rowIndexNearestXI][1 + index];
			double deltaPower = pow(delta, index);
			double fact = factorial(index);
			double multiplicant;
			double multiplicantContainer;
			for (int i = 0; i < index; i++) {
				multiplicant = x - differenceTable[rowIndexNearestXI - i][0];
				if (i != 0) {
					multiplicant = multiplicant*multiplicantContainer;
				}
				multiplicantContainer = multiplicant;
			}

			newValue = (diff / (deltaPower*fact))*multiplicantContainer;

			sumOfTerms = sumOfTerms + newValue;
			index++;
		}
	}
	else if (direction == DifferenceTableDirection::Center) {
		//for the central differences case, if there are even rows means there are two central rows 
		bool isEven = (intepolationTblRows % 2) == 0;

		//so find the x0 nearest to x from the last column of the DifferenceTable.
		int rowIndexNearestXI = 0;

		for (int i = 0; i < intepolationTblRows; i++) {
			if (differenceTable[i][differenceTblCols - 1] == 1) {
				rowIndexNearestXI = i;
				break;
			}
		}

		//if iseven and the selected row is actually not the second central row but a row higher in the table, then reduce the iteration count 
		//for going rightwards in the difference table for using higher order differences.
		int reduceTheCounter = 0;
		if (isEven) {
			int secondCentralRowIndex = intepolationTblRows / 2;
			if (rowIndexNearestXI < secondCentralRowIndex)
				reduceTheCounter = secondCentralRowIndex - rowIndexNearestXI;
		}

		int index = 1;
		double newValue;
		double x0 = differenceTable[rowIndexNearestXI][0];
		double x1 = differenceTable[rowIndexNearestXI + 1][0];
		double delta = x1 - x0;
		sumOfTerms = differenceTable[rowIndexNearestXI][1];
		//for the selected row x0, how far right in the row you can go inside the differencetable to use higher order differences.
		while (index < intepolationTblRows - rowIndexNearestXI - reduceTheCounter) {

			double diff = differenceTable[rowIndexNearestXI][1 + index];
			double deltaPower = pow(delta, index);
			double fact = factorial(index);
			double multiplicant;
			double multiplicantContainer;
			for (int i = 0; i < index; i++) {
				multiplicant = x - differenceTable[rowIndexNearestXI + i][0];
				if (i != 0) {
					multiplicant = multiplicant*multiplicantContainer;
				}
				multiplicantContainer = multiplicant;
			}

			newValue = (diff / (deltaPower*fact))*multiplicantContainer;

			sumOfTerms = sumOfTerms + newValue;
			index++;
		}
	}

	for (int i = 0; i < intepolationTblRows; i++)
		free(differenceTable[i]);
	free(differenceTable);

	return sumOfTerms;
}

//interpolation: when we have a value table of values of x and y, where y is f(x) but we do not know the f(x), then 
//to calculate a value of y for any other given x is called interpolation or extrapolation.
//	x		y
//	---------
//	1		4
//	2		5
//	4		11
//here, if we have to calculate the value of y for x=3 it is called interpolation and y for x=5 is called extrapolation.
//the word interpolate means estimate. 
//Estimation can be done by regression or curve fitting. In curve fitting problems, the constraint that the
//interpolant(estimate) has to go exactly through the data points is relaxed. It is only required to approach the 
//data points as closely as possible (within some other constraints). This requires parameterizing the potential 
//interpolants and having some way of measuring the error. In the simplest case this leads to least squares approximation.

//for the interpolation table given above, the spacing between x is not equal nor is the slope. Hence in this can't use LinearInterpolation
//(it needs equal slope) or NewtonInterpolation (as it needs equal spacing). For this case lagrange interpolation is used.
//explained at https://www3.nd.edu/~coast/jjwteach/www/www/30125/pdfnotes/lecture3_6v13.pdf. Also in git LagrangePolynomialInterpolation.pdf.
//try to see if the lagrange formula also be expressed in a somewhat 'Taylor series' form.
double LagrangePolynomialInterpolation(double x, double **interpolationTable, int intepolationTblRows) {
	throw NotImplementedException("Lagrange interpolation not implemented yet.");
}