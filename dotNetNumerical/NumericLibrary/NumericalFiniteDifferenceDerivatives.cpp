#include "stdafx.h"
#include "NumericLibrary.h"

//numerical derivatives using finite differences is used when the function f(x) iteself is not known. So we can't
//add a value delta to do something like (f(x0 + delta) - f(x0)) / ((x0 + delta) - x0). We have some values of x and f(x)
//known to us but not the function f(x)
//	x		y
//	---------
//	1		4
//	2		5
//	3		7
//	4		12
//	5		15
//so in this case we can't pick our small delta value to find the derivative but have to use the delta given to us 
//in the table. So if are calculating f(x) for 2.5, then we should pick x0=2(and use forward difference table) and we would have:
//f(x) = f(x0) + (f'(x0)/1!)(x- x0) + (f''(x0)/2!) (x- x0)^2 + (f'''(x0) / 3!) (x - x0) ^ 3
//now note here that x-x0=2.5-2=.5. Lets call it step size h which is the difference for value x for which we have
//to find f(x) from x0. this is usually different from delta used for finding derivatives
//now to  find f'(x0) using the values we have in the table, it would be:
//f'(x0)=(f(x0 + delta) - f(x0)) / ((x0 + delta) - x0). Now in this case we cant pick our own very small delta 
//for finding numerical derivatives and have to use spacing given to us in the table. We can't use arbitrary delta 
//as we do not know the function f(x). So we can write it as:
//(f(x0 + 1) - f(x0)) / ((x0 + 1) - x0) where 1 is the delta we have to use and the result of this can be found 
//from the difference table. Now more the points avaiable to us in the interpolation table, more the number 
//of higher order differences available to use. 
//

//Some questions arise here:
//1) Find out the mapping between the difference table so that derivative of any order can be found from it? Shown below
//2) What if the interpolation table does not equal spacing(for x values), then how do we find derivatives? We need equal 
//   spacing in the difference table for using it for interpolation as seen in the Newton method for polynomial interpolation.
//	 But do we need it for finding derivatives as well? I think yes as otherwise you would note that delta for the + and - part 
//	 in centered difference approximatin would be different. It would also make a difference for higher order differences in forward
//	 and backward difference table. A method called newtons divided difference should be used in case of x not being equi-spaced.
//	 I have not implemented it.
//3) This method is used when the function f(x) itself is not known. But if the function f(x) is known, in which case we case choose our own 
//	 very small delta value but how do we come to use more points than the oder of derivative. So for centered difference approximation, when 
//	 the function f(x) is know, the first order can derivative could be written as for better approximation:
//	 f'(x0) = -f(x0 + 2 * delta) + 8 * f(x0 + delta) - 8 * f(x0 - delta) + f(x0 - 2 * delta)) / (12 * delta)
//	 INSTEAD OF f'(x0)=(f(x0 + delta) - f(x0 - delta)) / ((x0 + delta) - (x0 - delta))
//	 so how did we come to use more than 2 points even for the first order derivative

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
//http://www.rsmas.miami.edu/personal/miskandarani/Courses/MSC321/lectfiniteDifference.pdf or in git DerivativesFromFiniteDifferences2.pdf

DifferenceTableDirection ComputeDifferenceTable(double **interpolationTable, double **differenceTable, int tableRows, int tableCols, double x, DifferenceTableDirection direction) {
	//for the central differences case, if there are even rows means there are two central rows 
	bool isEven = (tableRows % 2) == 0;

	//copy the 2 columns of interpolationTable to differenceTable
	for (int row = 0; row < tableRows; row++) {
		differenceTable[row][0] = interpolationTable[row][0];
		differenceTable[row][1] = interpolationTable[row][1];

		//encode distance from x of the current 'x' in table
		differenceTable[row][tableCols - 2] = abs(x - interpolationTable[row][0]);
	}


	if (direction == DifferenceTableDirection::Evaluate) {
		//find out what is x nearer to, top of the table or the bottom of the table or center
		double distanceFromTop = abs(x - interpolationTable[0][0]);
		double distanceFromBottom = abs(x - interpolationTable[tableRows - 1][0]);

		//if even rows, then two central rows whose x value has to compared against
		double distanceFromCenter;

		if (isEven) {
			double distanceFromFirstCentralRow = abs(x - interpolationTable[tableRows / 2][0]);
			double distanceFromSecondCentralRow = abs(x - interpolationTable[(tableRows / 2) + 1][0]);

			distanceFromCenter = distanceFromSecondCentralRow;
			if (distanceFromFirstCentralRow < distanceFromSecondCentralRow) {
				distanceFromCenter = distanceFromFirstCentralRow;
			}
		}
		else
			distanceFromCenter = abs(x - interpolationTable[(tableRows - 1) / 2][0]);

		if (distanceFromTop < distanceFromCenter && distanceFromTop < distanceFromBottom) {
			direction = DifferenceTableDirection::Forward;
		}
		else if (distanceFromBottom < distanceFromTop && distanceFromBottom < distanceFromCenter) {
			direction = DifferenceTableDirection::Backward;
		}
		else {
			direction = DifferenceTableDirection::Center;
		}
	}

	if (direction == DifferenceTableDirection::Forward) {
		//using the 'distance from x' value encoded above, find the index of nearest row to x and write 1 
		//against it. Set other rows values to 0
		double distance = DBL_MAX;
		int rowIndexNearestXI = 0;
		for (int i = 0; i < tableRows; i++) {
			if (differenceTable[i][tableCols - 2] < distance) {
				distance = differenceTable[i][tableCols - 2];
				differenceTable[rowIndexNearestXI][tableCols - 1] = 0;

				rowIndexNearestXI = i;
				differenceTable[rowIndexNearestXI][tableCols - 1] = 1;
			}
			else
				differenceTable[i][tableCols - 1] = 0;
		}

		for (int differenceOrder = 1; differenceOrder < tableRows; differenceOrder++) {
			for (int rowNo = 0; rowNo < tableRows - differenceOrder; rowNo++) {
				differenceTable[rowNo][1 + differenceOrder] = differenceTable[rowNo + 1][differenceOrder] - differenceTable[rowNo][differenceOrder];
			}
		}
	}
	else if (direction == DifferenceTableDirection::Backward) {
		//using the 'distance from x' value encoded above, find the index of nearest row to x and write 1 
		//against it. Set other rows values to 0
		double distance = DBL_MAX;
		int rowIndexNearestXI = tableRows - 1;
		for (int i = tableRows - 1; i >= 0; i--) {
			if (differenceTable[i][tableCols - 2] < distance) {
				distance = differenceTable[i][tableCols - 2];
				differenceTable[rowIndexNearestXI][tableCols - 1] = 0;

				rowIndexNearestXI = i;
				differenceTable[rowIndexNearestXI][tableCols - 1] = 1;
			}
			else
				differenceTable[i][tableCols - 1] = 0;
		}

		for (int differenceOrder = 1; differenceOrder < tableRows; differenceOrder++) {
			for (int rowNo = tableRows - 1; rowNo > differenceOrder - 1; rowNo--) {
				differenceTable[rowNo][1 + differenceOrder] = differenceTable[rowNo][differenceOrder] - differenceTable[rowNo - 1][differenceOrder];
			}
		}
	}
	else if (direction == DifferenceTableDirection::Center) {
		//using the 'distance from x' value encoded above, find the index of nearest row to x and write 1 
		//against it. Set other rows values to 0
		double distance = DBL_MAX;

		int upperHalfEnd = tableRows / 2;
		int lowerHalfStart = tableRows - upperHalfEnd;

		//in case of even number of rows, my central difference table creation logic always takes second row as selected row to populate the 
		//difference table. 
		//for (int differenceOrder = 1; differenceOrder < tableRows; differenceOrder++) {
		//	for (int rowNo = upperHalfEnd; rowNo < tableRows - differenceOrder; rowNo++) {
		//		differenceTable[rowNo][1 + differenceOrder] = differenceTable[rowNo + 1][differenceOrder] - differenceTable[rowNo - 1][differenceOrder];
		//		//how far have we gone down. use the down depth as a indicator for going up as well for calculating the differences
		//		int goUpCount = (rowNo - upperHalfEnd) * 2;
		//		if (goUpCount > 0) {
		//			differenceTable[rowNo - goUpCount][1 + differenceOrder] = differenceTable[rowNo - goUpCount + 1][differenceOrder] - differenceTable[rowNo - goUpCount - 1][differenceOrder];
		//		}
		//		if (isEven) //we are can go up one more row for cases when number of rows are even as there are two center rows
		//			differenceTable[rowNo - goUpCount - 1][1 + differenceOrder] = differenceTable[rowNo - goUpCount][differenceOrder] - differenceTable[rowNo - goUpCount - 2][differenceOrder];
		//	}
		//}

		for (int differenceOrder = 1; differenceOrder < tableRows - differenceOrder; differenceOrder++) {
			for (int i = differenceOrder; i < tableRows - differenceOrder; i++) {
				differenceTable[i][1 + differenceOrder] = differenceTable[i + 1][differenceOrder] -
															differenceTable[i - 1][differenceOrder];
			}
		}

		int rowIndexNearestXI = 0;

		for (int i = lowerHalfStart; i < tableRows; i++) {
			if (differenceTable[i][tableCols - 2] < distance) {
				distance = differenceTable[i][tableCols - 2];
				differenceTable[rowIndexNearestXI][tableCols - 1] = 0;

				rowIndexNearestXI = i;
				differenceTable[rowIndexNearestXI][tableCols - 1] = 1;
			}
			else
				differenceTable[i][tableCols - 1] = 0;
		}

		//if iseven, then that means both the ends are same. So reduce the upperend by -1
		isEven ? upperHalfEnd-- : upperHalfEnd;
		for (int i = upperHalfEnd; i >= 0; i--) {
			if (differenceTable[i][tableCols - 2] < distance) {
				distance = differenceTable[i][tableCols - 2];
				differenceTable[rowIndexNearestXI][tableCols - 1] = 0;

				rowIndexNearestXI = i;
				differenceTable[rowIndexNearestXI][tableCols - 1] = 1;
			}
			else
				differenceTable[i][tableCols - 1] = 0;
		}

	}

	return direction;
}
