// NumericalComputation.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"

std::map<double(*)(double), std::string> function_map;
#define REGISTER_FUNCTION(f) function_map[f] = BOOST_PP_STRINGIZE(f);

//4sinx-e^x=0
//solve for x=sin^-1(e^x/4)
//for iterative solver using successive approximation we have to rewrite the equation to get x on the LHS
double FourSinxMinusexRestructured(double x) {
	return asin(exp(x) / 4);
}

//4sinx-e^x=0
//for Newton Raphson method we DO NOT HAVE TO RESTRUCTURE the equation to x on the LHS.
//FourSinxMinusexRestructured solved by successive approximation iterative solver should give same result as FourSinxMinusex solved by Newton Raphson iterative solver
double FourSinxMinusex(double x) {
	return 4 * sin(x) - exp(x);
}

//x^3+12.1x^2+13.1x+22.2=0
//for Newton Raphson method we DO NOT HAVE TO RESTRUCTURE the equation to x on the LHS.
//FourSinxMinusexRestructured solved by successive approximation iterative solver should give same result as FourSinxMinusex solved by Newton Raphson iterative solver
double xCubeEquation(double x) {
	return pow(x, 3) + 12.1*pow(x, 2) + 13.1*x + 22.2;
}


double negativeExp(double x)
{
	return exp(-x);
}

//y=(1-x^2)^.5
double OneMinusXSq(double x)
{
	return pow((1 - pow(x, 2)), .5);
}

//y=1/(1+x)^.5
double OnePlusXSqrt(double x)
{
	return 1 / pow(1 + x, .5);
}

//polynomial f(x) = 3x ^ 3 + 2x ^ 2 + x - 2
double polynomialFunc(double x)
{
	return 3 * pow(x, 3) + 2 * pow(x, 2) + x - 2;
}

int main() {
	/*REGISTER_FUNCTION(negativeExp);
	printf("\nFunction name is: %s", function_map[negativeExp].c_str());*/

	//Get the machine precision
	GetMachinePrecision();
	printf("\n");

	//initial value is a guess and generally it is safe to use 0 as initial value
	//this is a iterative solver using successive approximation
	//more efficient iterative solvers also exist like the newton-raphson method
	double initialValue = 0;
	REGISTER_FUNCTION(FourSinxMinusex);
	IterativeRoot(FourSinxMinusex, initialValue);
	printf("\n");

	//calculate sin(x) where x =.25 using maclaurin series
	REGISTER_FUNCTION(sin);
	MaclaurinSeriesUsingNumericalDerivatives(sin, .25);
	sinxMaclaurinSeriesUsingDerivativesFromTigrometry(.25);
	printf("\n");

	MaclaurinSeriesUsingNumericalDerivatives(sin, .45);
	sinxMaclaurinSeriesUsingDerivativesFromTigrometry(.45);
	printf("\n");

	//now try maclaurin series for other functions

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	REGISTER_FUNCTION(negativeExp);
	MaclaurinSeriesUsingNumericalDerivatives(negativeExp, .2);
	printf("\n");

	//taylor series where x0=1, f(x) is log x and  we have to solve for f(x) where x is 2
	REGISTER_FUNCTION(log);
	TaylorSeriesUsingNumericalDerivatives(log, 1, 2);
	printf("\n");

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	REGISTER_FUNCTION(OneMinusXSq);
	MaclaurinSeriesUsingNumericalDerivatives(OneMinusXSq, .5);
	printf("\n");

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	REGISTER_FUNCTION(OnePlusXSqrt);
	MaclaurinSeriesUsingNumericalDerivatives(OnePlusXSqrt, 1 / pow(1.4, .5));
	printf("\n");

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	REGISTER_FUNCTION(cos);
	MaclaurinSeriesUsingNumericalDerivatives(cos, 1.5);
	printf("\n");

	//synthetic substitution can be used for solving polynomials. See the 'computational science' word doc. here we are passing the orderOfDerivative parameter set to 0. 
	//that means just calculate the polynomial.
	std::vector<double> polynomialCoefficients = { 3,2,1,-2 };
	double dsynsub0 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 0);
	printf("The solution of polynomial calculated using synthetic substitution is %lf\n", dsynsub0);
	dsynsub0 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 3, 0);
	printf("The solution of polynomial calculated using synthetic substitution is %lf\n", dsynsub0);

	//if the function is a polynomial, then instead of numerical derivaltives, we can also use Synthetic Substitution Derivatives
	//this catually gives us the derivative/factorial. So to get the exact derivative value, you need to, multiply with factorial value
	double dsynsub1 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 1) * factorial(1);
	double dsynsub2 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 2) * factorial(2);
	double dsynsub3 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 3) * factorial(3);

	//double dnum0 = derivative(polynomialFunc, 2, 0);
	double dnum1 = derivative(polynomialFunc, 2, 1);
	double dnum2 = derivative(polynomialFunc, 2, 2);
	double dnum3 = derivative(polynomialFunc, 2, 3);

	printf("The derivatives calculated using synthetic substitution are %lf, %lf and %lf\n", dsynsub1, dsynsub2, dsynsub3);
	printf("The derivatives calculated using numerical methods are %lf, %lf and %lf\n", dnum1, dnum2, dnum3);

	//taylor series where x0=2, f(x) is a polynomial f(x) = 3x^3 + 2x^2 + x - 2 and  we have to solve for f(x) where x is 3
	TaylorSeriesUsingSyntheticSubstitutionDerivatives(polynomialCoefficients, 2, 3);

	//taylor series where x0=2, f(x) is a polynomial f(x) = 3x^3 + 2x^2 + x - 2 and  we have to solve for f(x) where x is 3
	//polynomialFunc and polynomialCoefficients denote same polynomial
	REGISTER_FUNCTION(polynomialFunc);
	TaylorSeriesUsingNumericalDerivatives(polynomialFunc, 2, 3);
	printf("\n");

	//similarly calculate MaclaurinSeriesUsingNumericalDerivatives and MaclaurinSeriesUsingSyntheticSubstitution

	//difference table. This looks something that can be used to check stability . Not sure though. 
	//the calculation of difference table can be coded though.

	//for Newton Raphson method we DO NOT HAVE TO RESTRUCTURE the equation to x on the LHS.
	//FourSinxMinusexRestructured solved by successive approximation iterative solver should give same result as FourSinxMinusex solved by Newton Raphson iterative solver
	//note that the difference between Newton Raphson methods and successive approximation is that for successive approximation we restructure the function equation itself to get x on LHS
	//whereas for Newton Raphson method, we modify the taylor series for the function to get x on the LHS
	IterativeRoot(FourSinxMinusexRestructured, initialValue);
	IterativeRootNewtonRaphson(FourSinxMinusex, initialValue);
	printf("\n");

	//0 as the initial guess wont each time. So here for x^3+12.1x^2+13.1x+22.2=0, iterative root could not be found using initial value 0
	//but using -11, we found the root. There are ways to make educated guesses for initial values. page 299 of C.Xavier book on numerical methods explains that.
	IterativeRootNewtonRaphson(xCubeEquation, initialValue);
	IterativeRootNewtonRaphson(xCubeEquation, -11);
	initialValue = initialGuess(xCubeEquation);
	IterativeRootNewtonRaphson(xCubeEquation, initialValue);

	printf("\n");

	//Get some code timings
	cpu_timer timer;

	initialValue = initialGuess(xCubeEquation);
	IterativeRootNewtonRaphson(xCubeEquation, initialValue);

	timer.stop();

	cpu_times times = timer.elapsed();
	printf("Timings for evaluating IterativeRootNewtonRaphson %s\n", timer.format().c_str());
	printf("\n");

	//eigen library
	MatrixXd m(2, 2);
	m(0, 0) = 3;
	m(1, 0) = 2.5;
	m(0, 1) = -1;
	m(1, 1) = m(1, 0) + m(0, 1);
	printf("The determinant using Eigen library is %lf\n", m.determinant());
	std::cout << m << std::endl;
	printf("\n");

	//matrix manipulation
	//case 1
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

	printf("\n");

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

	times = timer.elapsed();
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

	/*printf("The input matrix for evaluating determinant is:\n");
	std::cout << matE << std::endl;*/

	timer.start();
	printf("The determinant using Eigen library is: %lf\n", matE.determinant());
	timer.stop();

	times = timer.elapsed();
	printf("Timings for evaluating determinant of the matrix using eigen library is %s\n", timer.format().c_str());
	printf("\n");

	//Determinant of the matrix using pivotal condensation method.
	matrixRows = 3;
	matrixColumns = 3;
	seed = 1;
	matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

	Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

	printf("The input matrix for evaluating determinant using pivotal condensation method is:\n");
	printMatrix(matA, matrixRows, matrixColumns);
	det = PivotalCondensationMethod(matA, matrixRows, matrixColumns);
	printf("And the determinant of the matrix using pivotal condensation method is: %lf\n", det);
	for (int i = 0; i < matrixRows; i++)
		free(matA[i]);
	free(matA);
	printf("\n");

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

	//timing comparisonsbetween pivotal condensation determinant evaluation and my naive determinant method
	matrixRows = 10;
	matrixColumns = 10;
	seed = 1;
	matA = Create2DimensionalMatrix(matrixRows, matrixColumns);

	Init2DimensionalMatrix(matA, matrixRows, matrixColumns, seed);

	//printf("The input matrix for evaluating determinant using pivotal condensation method is:\n");
	//printMatrix(matA, matrixRows, matrixColumns);

	timer.start();
	det = PivotalCondensationMethod(matA, matrixRows, matrixColumns);
	timer.stop();
	times = timer.elapsed();
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
	times = timer.elapsed();
	printf("Timings for evaluating determinant of the matrix using my code %s", timer.format().c_str());
	printf("And the determinant of the matrix using naive determinant method is: %lf\n", det);

	for (int i = 0; i < matrixRows; i++)
		free(matA[i]);
	free(matA);
	printf("\n");


	getchar();
	return 0;
}

