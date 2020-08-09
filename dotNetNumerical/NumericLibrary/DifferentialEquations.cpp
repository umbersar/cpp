#include "stdafx.h"
#include "NumericLibrary.h"

//Initial Value problem: The differential equation has to be first order differential equation. If it is not, then convert
//the higher order differential equation to a system of siumlataneous first order differential equations. And then solve
//the system of simulatenous equations using any of the initial value diff. eq. solvers
//fᵢ₊₁(x) = fᵢ(x) + f'ᵢ(x)*h. Given initial values of both x and y and the differential equation f'ᵢ(x), we can select a small h
//and keep iterating over the equation to plot the function as well as find the value of y at a particualar point

//taylor series is f(x) = f(x0) + (f'(x0)/1!)(x- x0) + (f''(x0)/2!) (x- x0)^2 + (f'''(0)/3!) (x- x0)^3 + ...
//truncating it to 2 terms we get the approximate series.
//f(x) ≈ f(x0) + (f'(x0)/1!)(x- x0)
//Now for eulers method, we are given the x0 and f(x0) and the derivative dy/dx equation(for eg. x+y). 
//We have to find y or f(x) for a particular x. So here we have to choose a small step h such that h=x-x0. Remember the difference 
//between step size h and delta. If x0=4 and the desired f(x)=23, then h!=23-4 but a small increment in x such as h=5-4.
//For a given h and x0, we can now calculate the f(x) by plugging the values of x and y at this point into the derivate equation.
//And we can keep repeating this step until we reach the desired f(x). Since we are iterating till we reach the desired x and
//its corresponding y, we can also write the equation in a iterative manner: fᵢ₊₁(x) = fᵢ(x) + f'ᵢ(x)*h

//The step size should be very small to get a more accurate approximation of the function.The smaller the step size, more would 
//be the number of steps required to get an answer at a given point or to draw the curve till a given point.This large number of 
//steps entails a high computational cost.For this reason, people usually employ alternative, higher-order approximation methods 
//such as Runge–Kutta methods or linear multistep methods, especially if a high accuracy is desired.The step size should be very small 
//to get a more accurate approximation of the function. Also, as the step size becomes smaller and smaller, you start inducing 
//error due to inability of the computer to represent them very small numbers. 
//How to improve eulers method ? To improve the solution we could either use smaller steps or use higher order approximation(not higher 
//order derivative) schemes which take into account the slopes at midway points as well for prediction so that for  large step, we are not way off the mark.

//eulers method and runge kutta method are initial value problem solvers involving differenetial equations.

double EulersMethod(double x0, double y0, double x, double(*f)(double, double)) {
	double h = sqrt(DBL_EPSILON);
	double diff = x - x0;

	// value of h has to be choosen to be small and should near perfectly be a multiple of diff
	int nSteps = abs(diff / h);
	//as explained here http://stackoverflow.com/q/39966190/364084, there is a bug in fmod.So instead do what its doco says it will do.
	//double rem = fmod(diff, h);// remainder(diff, h);
	double rem = diff - static_cast<int>(diff / h)*h;

	//rem already has the sign of diff. use that sign for h as well
	h = copysign(h, diff);

	double yi;

	for (int i = 0; i < nSteps; i++) {
		yi = y0 + f(x0, y0)*h;

		x0 = x0 + h;
		y0 = yi;
	}

	if (rem != 0.0) {
		yi = y0 + f(x0, y0)*rem;

		x0 = x0 + rem;
		if (abs(x0 - x) >= .0000001) {
			//should never get here
			throw NotComputableException("Solution for x=%d could not be computed. Change the step size to a more appropriate value.");
		}
		y0 = yi;
	}

	return yi;
}

double converge(double x0, double y0, double(*f)(double, double), double h) {
	//this is the predictor equation
	double oldY = y0 + f(x0, y0)  * h;
	double newY;

	//try to converge for a given number of times before giving up
	int index = 0;
	bool converged = false;
	while (index < 10) {
		//this is the corrector equation
		//recalculate yi using the mid point of slope as the actual slope value
		newY = y0 + (f(x0, y0) + f(x0 + h, oldY)) *(1.0 / 2) * h;

		if (abs(oldY - newY) <= .0001) {
			converged = true;
			break;
		}
		oldY = newY;
		index++;
	}

	if (!converged)
		throw NotComputableException("Eulers modified method did not converge.");

	return newY;
}

//eulers modified method is a type of predictor-corrector method
//even with a greater step size than eulers method, this gives better results.
//now now since the slope might change rapidly, use the of slope at current point to calculate a y value
//and compare that y value with the one calculated using average of slope calculated at current point and the next point 
//Keep doing it till the y value converges for the next point and then move ahead to the next step.
double EulersMethodModified(double x0, double y0, double x, double(*f)(double, double)) {
	double h = .1;
	double diff = x - x0;

	// value of h has to be choosen to be small and should near perfectly be a multiple of diff
	int nSteps = abs(diff / h);
	//as explained here http://stackoverflow.com/q/39966190/364084, there is a bug in fmod.So instead do what its doco says it will do.
	//double rem = fmodl(diff, h);//rem = remainderl(diff, h);
	double rem = diff - static_cast<int>(diff / h)*h;

	//rem already has the sign of diff. use that sign for h as well
	h = copysign(h, diff);

	double yi;

	for (int i = 0; i < nSteps; i++) {
		//x0 and y0 are initial values for each step in the iteration with yi being the solution at the end of a particular step/iteration.
		y0 = converge(x0, y0, f, h);
		x0 = x0 + h;
		yi = y0;
	}

	if (rem != 0.0) {
		yi = converge(x0, y0, f, rem);
	}

	return yi;
}

//runge kutta is more accurate and requires calculating the function only at some function points.
//this is the fourth order runge kutta where 4 referes to the order of aapproximation or accuracy
double RungeKuttaFourthOrder(double x0, double y0, double x, double(*f)(double, double)) {
	double h = .1;
	double diff = x - x0;

	// value of h has to be choosen to be small and should near perfectly be a multiple of diff
	int nSteps = abs(diff / h);
	//as explained here http://stackoverflow.com/q/39966190/364084, there is a bug in fmod.So instead do what its doco says it will do.
	//double rem = fmod(diff, h);// remainder(diff, h);
	double rem = diff - static_cast<int>(diff / h)*h;

	//rem already has the sign of diff. use that sign for h as well
	h = copysign(h, diff);

	double yi;
	for (int i = 0; i < nSteps; i++) {

		double k1 = h*f(x0, y0);
		double k2 = h*f(x0 + (1.0 / 2)*h, y0 + (1.0 / 2)*k1);
		double k3 = h*f(x0 + (1.0 / 2)*h, y0 + (1.0 / 2)*k2);
		double k4 = h*f(x0 + h, y0 + k3);

		double k = (1.0 / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

		yi = y0 + k;

		x0 = x0 + h;
		y0 = yi;
	}

	if (rem != 0.0) {

		double k1 = h*f(x0, y0);
		double k2 = h*f(x0 + (1.0 / 2)*h, y0 + (1.0 / 2)*k1);
		double k3 = h*f(x0 + (1.0 / 2)*h, y0 + (1.0 / 2)*k2);
		double k4 = h*f(x0 + h, y0 + k3);

		double k = (1.0 / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
		yi = y0 + k;

		x0 = x0 + rem;
		if (abs(x0 - x) >= .0000001) {
			//should never get here
			throw NotComputableException("Solution for x=%d could not be computed. Change the step size to a more appropriate value.");
		}
		y0 = yi;
	}
	return yi;
}

//I might make the runge kutta method more like a predictor corrector method where for each step, I compare the value by first 
//taking a particular h and then by taking h/2. If the result is not same I again use (h/2)/2 until we get similar values.
double RungeKuttaModified(double x0, double y0, double x, double(*f)(double, double)) {
	throw NotImplementedException("RungeKuttaModified not implemented as yet.");
}

//this should be the RungeKutta method but for keeping things easy for the reader of code, I have it in a separate method
//the fourth order here refers to 4rth order approximation, not 4rth order diff. eqs. To use any of the intial value diff. eq. solvers,
//the diff. eq has to be first order. If the diff. eq. is of a higher order, then convert the diff. eq. to a system of simulatenous first 
//order diff. eqs.
void RungeKuttaForTwoDifferentialEquationsFourthOrder(double *x, double *y, double *z, int sizeOfVariableArray, double(*f[])(double, double, double), double h) {
	double diff = x[sizeOfVariableArray - 1] - x[0];

	// value of h has to be choosen to be small and should near perfectly be a multiple of diff
	int nSteps = abs(diff / h);
	//as explained here http://stackoverflow.com/q/39966190/364084, there is a bug in fmod.So instead do what its doco says it will do.
	//double rem = fmod(diff, h);// remainder(diff, h);
	double rem = diff - static_cast<int>(diff / h)*h;

	//rem already has the sign of diff. use that sign for h as well
	h = copysign(h, diff);

	for (int i = 0; i < nSteps; i++) {

		double k1FirstEquation = h*f[0](x[i], y[i], z[i]);
		double k1SecondEquation = h*f[1](x[i], y[i], z[i]);

		double k2FirstEquation = h*f[0](x[i] + (1.0 / 2)*h, y[i] + (1.0 / 2)*k1FirstEquation, z[i] + (1.0 / 2)*k1SecondEquation);
		double k2SecondEquation = h*f[1](x[i] + (1.0 / 2)*h, y[i] + (1.0 / 2)*k1FirstEquation, z[i] + (1.0 / 2)*k1SecondEquation);

		double k3FirstEquation = h*f[0](x[i] + (1.0 / 2)*h, y[i] + (1.0 / 2)*k2FirstEquation, z[i] + (1.0 / 2)*k2SecondEquation);
		double k3SecondEquation = h*f[1](x[i] + (1.0 / 2)*h, y[i] + (1.0 / 2)*k2FirstEquation, z[i] + (1.0 / 2)*k2SecondEquation);

		double k4FirstEquation = h*f[0](x[i] + h, y[i] + k3FirstEquation, z[i] + k3SecondEquation);
		double k4SecondEquation = h*f[1](x[i] + h, y[i] + k3FirstEquation, z[i] + k3SecondEquation);

		double kFirstEquation = (1.0 / 6)*(k1FirstEquation + 2 * k2FirstEquation + 2 * k3FirstEquation + k4FirstEquation);
		double kSecondEquation = (1.0 / 6)*(k1SecondEquation + 2 * k2SecondEquation + 2 * k3SecondEquation + k4SecondEquation);

		//fill in the values of variables computed for the current step
		x[i + 1] = x[i] + h;
		y[i + 1] = y[i] + kFirstEquation;
		z[i + 1] = z[i] + kSecondEquation;
	}

	if (rem != 0.0) {
		double k1FirstEquation = rem*f[0](x[sizeOfVariableArray - 2], y[sizeOfVariableArray - 2], z[sizeOfVariableArray - 2]);
		double k1SecondEquation = rem*f[1](x[sizeOfVariableArray - 2], y[sizeOfVariableArray - 2], z[sizeOfVariableArray - 2]);

		double k2FirstEquation = rem*f[0](x[sizeOfVariableArray - 2] + (1.0 / 2)*h, y[sizeOfVariableArray - 2] + (1.0 / 2)*k1FirstEquation,
			z[sizeOfVariableArray - 2] + (1.0 / 2)*k1SecondEquation);
		double k2SecondEquation = rem*f[1](x[sizeOfVariableArray - 2] + (1.0 / 2)*h, y[sizeOfVariableArray - 2] + (1.0 / 2)*k1FirstEquation,
			z[sizeOfVariableArray - 2] + (1.0 / 2)*k1SecondEquation);

		double k3FirstEquation = rem*f[0](x[sizeOfVariableArray - 2] + (1.0 / 2)*h, y[sizeOfVariableArray - 2] + (1.0 / 2)*k2FirstEquation,
			z[sizeOfVariableArray - 2] + (1.0 / 2)*k2SecondEquation);
		double k3SecondEquation = rem*f[1](x[sizeOfVariableArray - 2] + (1.0 / 2)*h, y[sizeOfVariableArray - 2] + (1.0 / 2)*k2FirstEquation,
			z[sizeOfVariableArray - 2] + (1.0 / 2)*k2SecondEquation);

		double k4FirstEquation = rem*f[0](x[sizeOfVariableArray - 2] + h, y[sizeOfVariableArray - 2] + k3FirstEquation,
			z[sizeOfVariableArray - 2] + k3SecondEquation);
		double k4SecondEquation = rem*f[1](x[sizeOfVariableArray - 2] + h, y[sizeOfVariableArray - 2] + k3FirstEquation,
			z[sizeOfVariableArray - 2] + k3SecondEquation);

		double kFirstEquation = (1.0 / 6)*(k1FirstEquation + 2 * k2FirstEquation + 2 * k3FirstEquation + k4FirstEquation);
		double kSecondEquation = (1.0 / 6)*(k1SecondEquation + 2 * k2SecondEquation + 2 * k3SecondEquation + k4SecondEquation);

		//x[sizeOfVariableArray - 1] was set by the user. verify if adding the nsteps*h and rem  to starting x makes it equal to x[sizeOfVariableArray - 1]
		if (abs(x[sizeOfVariableArray - 1] - x[sizeOfVariableArray - 2] + rem) >= .0000001) {
			//should never get here
			throw NotComputableException("Solution could not be computed. Change the step size to a more appropriate value.");
		}

		//fill in the values of variables computed for the current step
		x[sizeOfVariableArray - 1] = x[sizeOfVariableArray - 2] + rem;
		y[sizeOfVariableArray - 1] = y[sizeOfVariableArray - 2] + kFirstEquation;
		z[sizeOfVariableArray - 1] = z[sizeOfVariableArray - 2] + kSecondEquation;
	}

}

//adamsbashforth is a predictor corrector method. We predict a value for a step h using predictor equation and then iterate 
//the corrector equatiion till we achieve the required tolerance.
//Let the interval over which we have to calculate the function values is x=x0 to xi. 
//for each small step h, Euler and rungekutta requires function info only at the begining of that step in the interval. 
//So we need (x0,y0) for y1 and (x1,y1) for y2.
//Whereas adamsbashforth requires function info not only at the begining of that step in the interval but also 3 prior
//function values. So we need (x(-3),y(-3)), (x(-2),y(-2)), (x0(-1),y(-1)) and (x0,y0) for y1.
//The predictor formula for this is the Fourth order 4-step explicit adams-bashforth method(note yᵢ is not fᵢ(x)):
//yᵢ₊₁ = yᵢ + (55*fᵢ(x) -59*fᵢ₋₁(x) + 37*fᵢ₋₂(x) - 9*fᵢ₋₃(x))*(h/24) 
//And the corrector formula is Fourth order 3 step Adams Moulton implicit method:
//yᵢ₊₁ = yᵢ + (9*fᵢ₊₁(x) + 19*fᵢ(x) - 5*fᵢ₋₁(x) + fᵢ₋₂(x))*(h/24)
//This method combines: Fourth order runge kutta(to calculate the prior 3 values) AND fourth order 4 step explicit adamsbasforth method
//And Fourth order 3 step Adams Moulton implicit method
//Fourth order refers to order of approximation(or accuracy). 
void AdamsBashforthMoultonFourthOrderForTwoDifferentialEquations(double x0, double y0, double z0, double x, double **xArr, double **yArr, double **zArr, int **sizeOfVariableArray, double(*f[])(double, double, double), double h) {
	int s = ceil((x - x0) / h) + 1 + 3;
	*sizeOfVariableArray = (int*)malloc(sizeof(int));
	**sizeOfVariableArray = s;
	*xArr = (double*)malloc(**sizeOfVariableArray * sizeof(double));
	*yArr = (double*)malloc(**sizeOfVariableArray * sizeof(double));
	*zArr = (double*)malloc(**sizeOfVariableArray * sizeof(double));
	xArr[0][3] = x0;
	yArr[0][3] = y0;
	zArr[0][3] = z0;
	xArr[0][**sizeOfVariableArray - 1] = x;


	double diff = xArr[0][**sizeOfVariableArray - 1] - xArr[0][3];

	// value of h has to be choosen to be small and should near perfectly be a multiple of diff
	int nSteps = abs(diff / h);
	//as explained here http://stackoverflow.com/q/39966190/364084, there is a bug in fmod.So instead do what its doco says it will do.
	//double rem = fmod(diff, h);// remainder(diff, h);
	double rem = diff - static_cast<int>(diff / h)*h;

	//rem already has the sign of diff. use that sign for h as well
	h = copysign(h, diff);

	//before we start adamsbashforth, we need to make sure we have the last 4 f values. If not, then use runge kutta(or Eulers) to find 
	//those values. For this reason runge kutta(or eulers) is called starter method. reverse the sign of h to go back
	int sizeOfPreviousValusesArrays = 4;
	double *xPreviousValArray = (double*)malloc(sizeOfPreviousValusesArrays * sizeof(double));
	double *yPreviousValArray = (double*)malloc(sizeOfPreviousValusesArrays * sizeof(double));
	double *zPreviousValArray = (double*)malloc(sizeOfPreviousValusesArrays * sizeof(double));
	//h = h*-1;//we are going in the other direction to compute prior values
	xPreviousValArray[0] = x0;
	yPreviousValArray[0] = y0;
	zPreviousValArray[0] = z0;
	xPreviousValArray[3] = -h * 4;
	RungeKuttaForTwoDifferentialEquationsFourthOrder(xPreviousValArray, yPreviousValArray, zPreviousValArray, sizeOfPreviousValusesArrays, f, h);

	//reverse the sign back to original
	//h = h*-1;

	//copy the data from PreviousValArrays in reverse. Use it to calculate the adams-bashforth logic. 
	for (size_t i = 0; i < 4; i++) {
		xArr[0][i] = xPreviousValArray[3 - i];
		yArr[0][i] = yPreviousValArray[3 - i];
		zArr[0][i] = zPreviousValArray[3 - i];
	}

	for (int i = 3; i < **sizeOfVariableArray - 1; i++) {
		//predictor formula for this is:  fᵢ₊₁(x) = fᵢ(x) + (55 * fᵢ(x) - 59 * fᵢ₋₁(x) + 37 * fᵢ₋₂(x) - 9 * fᵢ₋₃(x))*(h / 24)

		yArr[0][i + 1] = yArr[0][i] + (h / 24)*(55 * f[0](xArr[0][i], yArr[0][i], zArr[0][i]) - 59 * f[0](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1])
			+ 37 * f[0](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]) - 9 * f[0](xArr[0][i - 3], yArr[0][i - 3], zArr[0][i - 3]));

		zArr[0][i + 1] = zArr[0][i] + (h / 24)*(55 * f[1](xArr[0][i], yArr[0][i], zArr[0][i]) - 59 * f[1](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1])
			+ 37 * f[1](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]) - 9 * f[1](xArr[0][i - 3], yArr[0][i - 3], zArr[0][i - 3]));

		//corrector formula is: fᵢ₊₁(x) = fᵢ(x) + (9*fᵢ₊₁(x) + 19*fᵢ(x) - 5*fᵢ₋₁(x) + fᵢ₋₂(x))*(h/24)
		int index = 0;
		double yNewVal;
		double zNewVal;
		while (index < 10) {
			yNewVal = yArr[0][i] + (h / 24)*(9 * f[0](xArr[0][i + 1], yArr[0][i + 1], zArr[0][i + 1]) + 19 * f[0](xArr[0][i], yArr[0][i], zArr[0][i])
				- 5 * f[0](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1]) + f[0](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]));

			zNewVal = zArr[0][i] + (h / 24)*(9 * f[1](xArr[0][i + 1], yArr[0][i + 1], zArr[0][i + 1]) + 19 * f[1](xArr[0][i], yArr[0][i], zArr[0][i])
				- 5 * f[1](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1]) + f[1](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]));


			//check for tolerance
			if (abs(yNewVal - yArr[0][i + 1]) < .0001 && abs(zNewVal - zArr[0][i + 1]) < .0001) {
				yArr[0][i + 1] = yNewVal;
				zArr[0][i + 1] = zNewVal;
				break;
			}

			yArr[0][i + 1] = yNewVal;
			zArr[0][i + 1] = zNewVal;
			index++;
		}

		//fill in the value of variable for the current step
		xArr[0][i + 1] = xArr[0][i] + h;
	}

	if (rem != 0.0) {
		//predictor formula for this is:  fᵢ₊₁(x) = fᵢ(x) + (55 * fᵢ(x) - 59 * fᵢ₋₁(x) + 37 * fᵢ₋₂(x) - 9 * fᵢ₋₃(x))*(h / 24)
		int i = **sizeOfVariableArray - 1;

		yArr[0][i + 1] = yArr[0][i] + (rem / 24)*(55 * f[0](xArr[0][i], yArr[0][i], zArr[0][i]) - 59 * f[0](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1])
			+ 37 * f[0](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]) - 9 * f[0](xArr[0][i - 3], yArr[0][i - 3], zArr[0][i - 3]));

		zArr[0][i + 1] = zArr[0][i] + (rem / 24)*(55 * f[1](xArr[0][i], yArr[0][i], zArr[0][i]) - 59 * f[1](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1])
			+ 37 * f[1](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]) - 9 * f[1](xArr[0][i - 3], yArr[0][i - 3], zArr[0][i - 3]));

		//corrector formula is: fᵢ₊₁(x) = fᵢ(x) + (9*fᵢ₊₁(x) + 19*fᵢ(x) - 5*fᵢ₋₁(x) + fᵢ₋₂(x))*(h/24)
		int index = 0;
		double yNewVal;
		double zNewVal;
		while (index < 10) {
			yNewVal = yArr[0][i] + (rem / 24)*(9 * f[0](xArr[0][i + 1], yArr[0][i + 1], zArr[0][i + 1]) + 19 * f[0](xArr[0][i], yArr[0][i], zArr[0][i])
				- 5 * f[0](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1]) + f[0](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]));

			zNewVal = zArr[0][i] + (rem / 24)*(9 * f[1](xArr[0][i + 1], yArr[0][i + 1], zArr[0][i + 1]) + 19 * f[1](xArr[0][i], yArr[0][i], zArr[0][i])
				- 5 * f[1](xArr[0][i - 1], yArr[0][i - 1], zArr[0][i - 1]) + f[1](xArr[0][i - 2], yArr[0][i - 2], zArr[0][i - 2]));


			//check for tolerance
			if (yNewVal - yArr[0][i + 1] < .0001 && zNewVal - zArr[0][i + 1] < .0001) {
				yArr[0][i + 1] = yNewVal;
				zArr[0][i + 1] = zNewVal;
				break;
			}

			yArr[0][i + 1] = yNewVal;
			zArr[0][i + 1] = zNewVal;
			index++;
		}

		//x[sizeOfVariableArray - 1] was set by the user. verify if adding the nsteps*h and rem  to starting x makes it equal to x[sizeOfVariableArray - 1]
		if (abs(xArr[**sizeOfVariableArray - 1] - xArr[**sizeOfVariableArray - 2] + rem) >= .0000001) {
			//should never get here
			throw NotComputableException("Solution could not be computed. Change the step size to a more appropriate value.");
		}

		//fill in the value of variable for the current step
		xArr[0][i + 1] = xArr[0][i] + rem;
	}

}