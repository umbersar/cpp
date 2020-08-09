#include "stdafx.h"
#include "NumericLibrary.h"

//http://math.nyu.edu/~atm262/fall06/compmethods/a1/DerivativesApproximationByFiniteDifferences.pdf
//numerical derivative using centered difference approximation. derivative2 down below is also using 
//higher order centered difference approximation. Higher order approximation is different from higher order derivative
NUMERICLIBRARY_API double derivative1(double(*f)(double), double x0, int orderOfDerivative) {
	//const double delta = 1.0e-6;
	double delta = nextafter(x0, 1);
	//delta = sqrt(DBL_EPSILON);
	delta = 0.0625;//this delta also gave reasonable results

	if (orderOfDerivative == 1) {
		if (((x0 + delta) - (x0 - delta)) != 0.0) {
			double term1 = f(x0 + delta);
			double term2 = f(x0 - delta);

			if (isnan(term1) || isnan(term2)|| isinf(term1) || isinf(term2))
				throw NotComputableException("Floating point division by zero error caused inside the function whose derivative is being computed.");
			//return (f(x0 + delta) - f(x0 - delta)) / ((x0 + delta) - (x0 - delta));
			return (term1 - term2) / ((x0 + delta) - (x0 - delta));
		}
		else
			throw NotComputableException("Floating point division by zero error caused while calculating the derivative.");
	}
	else {
		if (((x0 + delta) - (x0 - delta)) != 0.0)
			return (derivative1(f, x0 + delta, orderOfDerivative - 1) - derivative1(f, x0 - delta, orderOfDerivative - 1))
			/ ((x0 + delta) - (x0 - delta));
		else
			throw NotComputableException("Floating point division by zero error caused while calculating the derivative.");
	}
}

//http://math.nyu.edu/~atm262/fall06/compmethods/a1/DerivativesApproximationByFiniteDifferences.pdf or in git DerivativesApproximationByFiniteDifferences.pdf
//this numerical derivate uses a higher order centered difference approximation. Note the difference bw higher order centered difference approximation  
//and higher order derivatives. By using more points, we are trying to aim for a higher order approximation(not higher order derivative).
//Question: For centered difference approximation, INSTEAD OF f'(x0)=(f(x0 + delta) - f(x0 - delta)) / ((x0 + delta) - (x0 - delta)),
//how did we derive the current formula? DerivativesFromFiniteDifferences2.pdf shows how we can derive this formula by using the taylor series.
//Todo: the delta of choice has to be found for different functions we are finding derivatives for: http://stackoverflow.com/a/637969/364084
//and http://stackoverflow.com/a/632538/364084. But I am hardcoding it to .0625 here.
NUMERICLIBRARY_API double derivative2(double(*f)(double), double x0, int orderOfDerivative) {
	//const double delta = 1.0e-6;
	double delta = nextafter(x0, 1);
	//delta = sqrt(DBL_EPSILON);
	delta = 0.0625;//this delta gave reasonable results

	if (orderOfDerivative == 1) {
		double term1 = f(x0 + 2 * delta);
		double term2 = f(x0 + delta);
		double term3 = f(x0 - delta);
		double term4 = f(x0 - 2 * delta);

		if (isnan(term1) || isnan(term2) || isnan(term3) || isnan(term4)
			|| isinf(term1) || isinf(term2) || isinf(term3) || isinf(term4))
			throw NotComputableException("Floating point division by zero error caused inside the function whose derivative is being computed.");
		//return (-f(x0 + 2 * delta) + 8 * f(x0 + delta) - 8 * f(x0 - delta) + f(x0 - 2 * delta)) / (12 * delta);
		return (-term1 + 8 * term2 - 8 * term3 + term4) / (12 * delta);
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
NUMERICLIBRARY_API double derivative(double(*f)(double), double x0, int orderOfDerivative) {
	return derivative2(f, x0, orderOfDerivative);
	//should I use epsilon here or nextafter value. For some reason, sqrt(DBL_EPSILON) provided 
	//the correct approximate result. 
	double delta = nextafter(x0, 1);
	//double delta = sqrt(DBL_EPSILON);
	delta = 0.0625;//this delta also gave reasonable results

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
NUMERICLIBRARY_API double derivativeDividedByFactorialUsingSyntheticSubstitution(std::vector<double> coefficients, double x, int orderOfDerivative) {
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

NUMERICLIBRARY_API double sinDerivativeFromTrignometry(double x0, int orderOfDerivative) {
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