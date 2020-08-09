#include "stdafx.h"

//http://math.nyu.edu/~atm262/fall06/compmethods/a1/DerivativesApproximationByFiniteDifferences.pdf
//numerical derivative using centered difference approximation. I have implemented it in a incorrectly.
double derivative1(double(*f)(double), double x0, int orderOfDerivative) {
	//const double delta = 1.0e-6;
	double delta = nextafter(x0, 1);
	delta = sqrt(DBL_EPSILON);

	if (orderOfDerivative == 1) {
		return (f(x0 + delta) - f(x0 - delta)) / ((x0 + delta) - (x0 - delta));
	}
	else {
		return (derivative1(f, x0 + delta, orderOfDerivative - 1) - derivative1(f, x0 - delta, orderOfDerivative - 1))
			/ ((x0 + delta) - (x0 - delta));
	}
}

//http://math.nyu.edu/~atm262/fall06/compmethods/a1/DerivativesApproximationByFiniteDifferences.pdf
//this numerical derivate using centered difference approximation implementation uses more points for calculating a derivative so that 
//we can use this implementation for higher order derivatives. The mathnetnumerics 
//code makes you use more points than the order of derivative. So I am replicating that behaviour.
double derivative2(double(*f)(double), double x0, int orderOfDerivative) {
	//const double delta = 1.0e-6;
	double delta = nextafter(x0, 1);
	delta = sqrt(DBL_EPSILON);//this delta gave reasonable results
	delta = 0.0625;//this delta gave reasonable results

	if (orderOfDerivative == 1) {
		return (-f(x0 + 2 * delta) + 8 * f(x0 + delta) - 8 * f(x0 - delta) + f(x0 - 2 * delta)) / (12 * delta);
	}
	else {
		return (-derivative2(f, x0 + 2 * delta, orderOfDerivative - 1) + 8 * derivative2(f, x0 + delta, orderOfDerivative - 1) - 8 * derivative2(f, x0 - delta, orderOfDerivative - 1) + derivative2(f, x0 - 2 * delta, orderOfDerivative - 1))
			/ (12 * delta);
	}
}

//http://math.nyu.edu/~atm262/fall06/compmethods/a1/DerivativesApproximationByFiniteDifferences.pdf
//this numerical derivative approximation is using finite difference method, specifically
//the forward difference approximation. mathnumerics, derivate1 and derivate2 use the centered
//difference method
double derivative(double(*f)(double), double x0, int orderOfDerivative) {
	//should I use epsilon here or nextafter value. For some reason, sqrt(DBL_EPSILON) provided 
	//the correct approximate result. 
	//double delta = nextafter(x0, 1);
	//double delta = sqrt(DBL_EPSILON);
	double delta = 0.0625;//this delta also gave reasonable results

	if (orderOfDerivative == 1) {
		return (f(x0 + delta) - f(x0)) / ((x0 + delta) - x0);
	}
	else {
		return (derivative(f, x0 + delta, orderOfDerivative - 1) - derivative(f, x0, orderOfDerivative - 1))
			/ (x0 + delta - x0);
	}
}

//for polynomials, we can also use synthetic substitution to calculate derivatives/factorial.
//when the order of the derivative is 0, this method gives a result for a polynomial using synthetic substitution 
//the result of this method needs to be multiplied by facorial to get the derivative value.
double derivativeDividedByFactorialUsingSyntheticSubstitution(std::vector<double> coefficients, double x, int orderOfDerivative) {
	//int array_size = sizeof(coefficients) / sizeof(double);

	std::vector<double> resultCoefficients;
	double result;
	for (size_t order = 0; order <= orderOfDerivative; order++) {
		result = 0;
		//resultCoefficients.push_back(0);
		for (size_t i = 0; i < coefficients.size(); i++) {
			result = (coefficients[i] + result);
			if (resultCoefficients.size() < coefficients.size() - 1)
				resultCoefficients.push_back(result);

			if (i < coefficients.size() - 1) {
				result = x * result;
			}

		}
		coefficients = resultCoefficients;
		resultCoefficients.clear();
	}
	return result;
}

double sinDerivativeFromTrignometry(double x0, int orderOfDerivative) {
	if (orderOfDerivative > 5) {
		printf("Trignometric derivates supported only till 5th order.\n");
		return NAN;
	}

	double(*f[10])(double);
	f[0] = cos;
	f[1] = sin;
	f[2] = cos;
	f[3] = sin;
	f[4] = cos;

	double derivativeValue = f[orderOfDerivative - 1](x0);

	if (orderOfDerivative == 2 || orderOfDerivative == 3) {
		derivativeValue = -derivativeValue;
	}

	return derivativeValue;
}