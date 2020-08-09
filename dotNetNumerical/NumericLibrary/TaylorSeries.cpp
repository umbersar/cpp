#include "stdafx.h"
#include "NumericLibrary.h"

NUMERICLIBRARY_API double MaclaurinSeriesUsingNumericalDerivatives(double(*f)(double), double x, double accuracy, int maxIterations) {
	//the tolerance or precision or accuracy to which you want to calculate the sin(x)

	int index = 1;

	double newValue;
	double lastValue = f(0);
	double sumOfTerms = lastValue;
	while (true) {
		if (index == maxIterations) {
			//todo: throw exception here
			printf("Could not solve to the required accuracy %lf in %d iterations. Relax the max iteration count to see if a solution is possible.\n", accuracy, index);
			printf("The result computed using numerical derivatives in %d iterations is %lf\n", index, sumOfTerms);
			return NAN;
		}

		//there can be 2 issues here. either the derivative function is wrong in which
		//case i can look for alternative implementations OR the delta used for numerical derivative 
		//calculation is not appropriate one. Remeber derivative is slope. We have to add a very small number
		//to x to calcualte the slope at that very point and that is delta.
		double der = derivative(f, 0, index);
		double xpower = pow(x, index);
		double fact = factorial(index);
		newValue = der * xpower / fact;

		sumOfTerms = sumOfTerms + newValue;

		if (fabs(newValue - lastValue) <= accuracy) {
			printf("Solution of using numerical derivatives for x=%lf is: %lf\n", x, sumOfTerms);
			return sumOfTerms;
		}
		else {
			lastValue = newValue;
		}

		index++;
	}
}

NUMERICLIBRARY_API double sinxMaclaurinSeriesUsingDerivativesFromTigrometry(double x, double accuracy, int maxIterations) {
	//the tolerance or precision to which you want to calculate the sin(x)

	int index = 1;

	double newValue;
	double lastValue = sin(0);
	double sumOfTerms = lastValue;
	while (true) {
		if (index == maxIterations) {
			//todo: throw exception here
			printf("Could not solve to the required accuracy %lf in %d iterations. Relax the max iteration count to see if a solution is possible.\n", accuracy, index);
			printf("The result computed using trignometric derivatives in %d iterations is %lf\n", index, sumOfTerms);
			return NAN;
		}

		double der = sinDerivativeFromTrignometry(0, index);
		double power = pow(x, index);
		double fact = factorial(index);
		newValue = der * power / fact;

		sumOfTerms = sumOfTerms + newValue;

		if (fabs(newValue - lastValue) <= accuracy) {
			printf("Solution of using trignometric derivatives for x=%lf is: %lf\n", x, sumOfTerms);
			return sumOfTerms;
		}
		else {
			lastValue = newValue;
		}

		index++;
	}
}

//solve f(x) =  sin(x) using maclaurin series
//Maclaurin series is one in which x0 is 0
//solve for f(x) = f(0) + (f'(0)/1!)(x) + (f''(0)/2!)(x)^2 + (f'''(0)/3!)(x)^3
NUMERICLIBRARY_API void sinxMaclaurin(double x) {
	//now we have to use maclaurin series(form of taylor series) to calculate sin(x).
	//now since we know that f'(x)=cos x and f''(x)=-sin x and f'''(x)=-cos x and so on.
	//so we do not have to calculate the derivatives. So for this excercise, we will solve this sin(x)
	//by first solving for derivatives ourselves and once by using the information that we have about derivative of sin(x).

	MaclaurinSeriesUsingNumericalDerivatives(sin, x);
	sinxMaclaurinSeriesUsingDerivativesFromTigrometry(x);
}

NUMERICLIBRARY_API double TaylorSeriesUsingNumericalDerivatives(double(*f)(double), double x0, double x, double accuracy, int maxIterations) {
	//the tolerance or precision or accuracy to which you want to calculate the sin(x)

	int index = 1;

	double newValue;
	double lastValue = f(x0);
	double sumOfTerms = lastValue;
	while (true) {
		if (index == maxIterations) {
			//todo: throw exception here
			printf("Could not solve to the required accuracy %lf in %d iterations. Relax the max iteration count to see if a solution is possible.\n", accuracy, index);
			printf("The result computed using numerical derivatives in %d iterations is %lf\n", index, sumOfTerms);
			return NAN;
		}

		//there can be 2 issues here. either the derivative function is wrong in which
		//case i can look for alternative implementations OR the delta used for numerical derivative 
		//calculation is not appropriate one. Remeber derivative is slope. We have to add a very small number
		//to x to calcualte the slope at that very point and that is delta.
		double der = derivative(f, x0, index);
		double xpower = pow(x - x0, index);
		double fact = factorial(index);
		newValue = der * xpower / fact;

		sumOfTerms = sumOfTerms + newValue;

		if (fabs(newValue - lastValue) <= accuracy) {
			printf("Solution of using numerical derivatives for x=%lf is: %lf\n", x, sumOfTerms);
			return sumOfTerms;
		}
		else {
			lastValue = newValue;
		}

		index++;
	}
}

NUMERICLIBRARY_API double TaylorSeriesUsingSyntheticSubstitutionDerivatives(std::vector<double> polynomialCoefficients, double x0, double x, double accuracy, int maxIterations) {
	//the tolerance or precision to which you want to calculate the sin(x)

	int index = 1;

	double newValue;
	double lastValue = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, x0, 0);
	double sumOfTerms = lastValue;
	while (true) {
		if (index == maxIterations) {
			//todo: throw exception here
			printf("Could not solve polynomial using synthetic substitution derivation to the required accuracy %lf in %d iterations. Relax the max iteration count to see if a solution is possible.\n", accuracy, index);
			printf("The result computed using numerical derivatives in %d iterations is %lf\n", index, sumOfTerms);
			return NAN;
		}

		//there can be 2 issues here. either the derivative function is wrong in which
		//case i can look for alternative implementations OR the delta used for numerical derivative 
		//calculation is not appropriate one. Remeber derivative is slope. We have to add a very small number
		//to x to calcualte the slope at that very point and that is delta.
		double der = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, x0, index);
		double xpower = pow(x - x0, index);
		double fact = factorial(index);
		newValue = der * xpower;// we do not divide by fact in case of synthetic substitution.

		sumOfTerms = sumOfTerms + newValue;

		if (fabs(newValue - lastValue) <= accuracy) {
			printf("\nSolution of polynomial using synthetic substitution derivation for x=%lf is: %lf", x, sumOfTerms);
			return sumOfTerms;
		}
		else {
			lastValue = newValue;
		}

		index++;
	}
}

//taylor series is f(x) = f(x0) + (f'(x0)/1!)(x- x0) + (f''(x0)/2!) (x- x0)^2 + (f'''(0)/3!) (x- x0)^3
//where x0 is the starting point from which we have to proceed ahead.
//now if x0 is 0 then it is Maclaurin series. Note that x-x0 is the step size and if we consider
//only the first two terms of the taylor series, then first term is the starting point on the curve and (x-x0)
// is the step size which is being multiplied by the slope at the point to get the next point on the curve
NUMERICLIBRARY_API void TaylorSeries(double(*f)(double), double x0, double x) {
	TaylorSeriesUsingNumericalDerivatives(f, x0, x);
}

//solve f(x) using maclaurin series
//Maclaurin series is taylor series in which x0 is 0
//solve for f(x) = f(0) + (f'(0)/1!)(x) + (f''(0)/2!)(x)^2 + (f'''(0)/3!)(x)^3
void MaclaurinSeries(double x) {
	MaclaurinSeriesUsingNumericalDerivatives(sin, x);
}