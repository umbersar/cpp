// NumericLibraryConsumer.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "..\NumericLibrary\NumericLibrary.h"

//for GNUPLOT
#include "gnuplot_i.hpp"
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <conio.h>   //for getch(), needed in wait_for_key()
#include <windows.h> //for Sleep()
void sleep(int i)
{
	Sleep(i * 1000);
}
#endif

#define SLEEP_LGTH 2  // sleep time in seconds
#define NPOINTS    50 // length of array

void wait_for_key(); // Programm halts until keypress
void GnuPlotSamples();


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

	GnuPlotSamples();
	/*REGISTER_FUNCTION(negativeExp);
	printf("\nFunction name is: %s", function_map[negativeExp].c_str());*/

	//Get the machine precision
	GetMachinePrecision();
	printf("\n");

	//initial value is a guess and generally it is safe to use 0 as initial value
	//this is a iterative solver using successive approximation
	//more efficient iterative solvers also exist like the newton-raphson method
	double initialValue = 0;
	IterativeRoot(FourSinxMinusexRestructured, initialValue);
	printf("\n");

	//calculate sin(x) where x =.25 using maclaurin series
	MaclaurinSeriesUsingNumericalDerivatives(sin, .25);
	sinxMaclaurinSeriesUsingDerivativesFromTigrometry(.25);
	printf("\n");

	MaclaurinSeriesUsingNumericalDerivatives(sin, .45);
	sinxMaclaurinSeriesUsingDerivativesFromTigrometry(.45);
	printf("\n");

	//now try maclaurin series for other functions

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	MaclaurinSeriesUsingNumericalDerivatives(negativeExp, .2);
	printf("\n");

	//taylor series where x0=1, f(x) is log x and  we have to solve for f(x) where x is 2
	TaylorSeriesUsingNumericalDerivatives(log, 1, 2);
	printf("\n");

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	MaclaurinSeriesUsingNumericalDerivatives(OneMinusXSq, .5);
	printf("\n");

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
	MaclaurinSeriesUsingNumericalDerivatives(OnePlusXSqrt, 1 / pow(1.4, .5));
	printf("\n");

	//todo: i get correct result when i set the "Solution Platforms" to mixed platforms but start getting wrong results 
	//when i set to x64
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

	//Gauss elimination method
	matrixRows = 3;
	matrixColumns = 3;
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

	getchar();
	return 0;
}

void wait_for_key()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
	std::cout << std::endl << "Press any key to continue..." << std::endl;

	FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
	_getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	std::cout << std::endl << "Press ENTER to continue..." << std::endl;

	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail());
	std::cin.get();
#endif
	return;
}


void GnuPlotSamples()
{
	//
	// Using the GnuplotException class
	//
	try
	{
		Gnuplot::set_GNUPlotPath("C:/gnuplot/bin");
		Gnuplot g1("lines");

		//
		// Slopes
		//
		std::cout << "*** plotting slopes" << std::endl;
		g1.set_title("Slopes\\nNew Line");

		std::cout << "y = x" << std::endl;
		g1.plot_slope(1.0, 0.0, "y=x");

		std::cout << "y = 2*x" << std::endl;
		g1.plot_slope(2.0, 0.0, "y=2x");

		std::cout << "y = -x" << std::endl;
		g1.plot_slope(-1.0, 0.0, "y=-x");
		g1.unset_title();

		//
		// Equations
		//
		g1.reset_plot();
		std::cout << std::endl << std::endl << "*** various equations" << std::endl;

		std::cout << "y = sin(x)" << std::endl;
		g1.plot_equation("sin(x)", "sine");

		std::cout << "y = log(x)" << std::endl;
		g1.plot_equation("log(x)", "logarithm");

		std::cout << "y = sin(x) * cos(2*x)" << std::endl;
		g1.plot_equation("sin(x)*cos(2*x)", "sine product");

		//
		// Styles
		//
		g1.reset_plot();
		std::cout << std::endl << std::endl << "*** showing styles" << std::endl;

		std::cout << "sine in points" << std::endl;
		g1.set_pointsize(0.8).set_style("points");
		g1.plot_equation("sin(x)", "points");

		std::cout << "sine in impulses" << std::endl;
		g1.set_style("impulses");
		g1.plot_equation("sin(x)", "impulses");

		std::cout << "sine in steps" << std::endl;
		g1.set_style("steps");
		g1.plot_equation("sin(x)", "steps");

		//
		// Save to ps
		//
		g1.reset_all();
		std::cout << std::endl << std::endl << "*** save to ps " << std::endl;

		std::cout << "y = sin(x) saved to test_output.ps in working directory" << std::endl;
		//      g1.savetops("test_output");
		g1.savetofigure("test_output.ps", "postscript color");
		g1.set_style("lines").set_samples(300).set_xrange(0, 5);
		g1.plot_equation("sin(12*x)*exp(-x)").plot_equation("exp(-x)");

		g1.showonscreen(); // window output


						   //
						   // User defined 1d, 2d and 3d point sets
						   //
		std::vector<double> x, y, y2, dy, z;

		for (unsigned int i = 0; i < NPOINTS; i++)  // fill double arrays x, y, z
		{
			x.push_back((double)i);             // x[i] = i
			y.push_back((double)i * (double)i); // y[i] = i^2
			z.push_back(x[i] * y[i]);           // z[i] = x[i]*y[i] = i^3
			dy.push_back((double)i * (double)i / (double)10); // dy[i] = i^2 / 10
		}
		y2.push_back(0.00);
		y2.push_back(0.78);
		y2.push_back(0.97);
		y2.push_back(0.43);
		y2.push_back(-0.44);
		y2.push_back(-0.98);
		y2.push_back(-0.77);
		y2.push_back(0.02);


		g1.reset_all();
		std::cout << std::endl << std::endl << "*** user-defined lists of doubles" << std::endl;
		g1.set_style("impulses").plot_x(y, "user-defined doubles");

		g1.reset_plot();
		std::cout << std::endl << std::endl << "*** user-defined lists of points (x,y)" << std::endl;
		g1.set_grid();
		g1.set_style("points").plot_xy(x, y, "user-defined points 2d");

		g1.reset_plot();
		std::cout << std::endl << std::endl << "*** user-defined lists of points (x,y,z)" << std::endl;
		g1.unset_grid();
		g1.plot_xyz(x, y, z, "user-defined points 3d");

		g1.reset_plot();
		std::cout << std::endl << std::endl << "*** user-defined lists of points (x,y,dy)" << std::endl;
		g1.plot_xy_err(x, y, dy, "user-defined points 2d with errorbars");


		//
		// Multiple output screens
		//
		std::cout << std::endl << std::endl;
		std::cout << "*** multiple output windows" << std::endl;

		g1.reset_plot();
		g1.set_style("lines");
		std::cout << "window 1: sin(x)" << std::endl;
		g1.set_grid().set_samples(600).set_xrange(0, 300);
		g1.plot_equation("sin(x)+sin(x*1.1)");

		g1.set_xautoscale().replot();

		Gnuplot g2;
		std::cout << "window 2: user defined points" << std::endl;
		g2.plot_x(y2, "points");
		g2.set_smooth().plot_x(y2, "cspline");
		g2.set_smooth("bezier").plot_x(y2, "bezier");
		g2.unset_smooth();

		Gnuplot g3("lines");
		std::cout << "window 3: log(x)/x" << std::endl;
		g3.set_grid();
		g3.plot_equation("log(x)/x", "log(x)/x");

		Gnuplot g4("lines");
		std::cout << "window 4: splot x*x+y*y" << std::endl;
		g4.set_zrange(0, 100);
		g4.set_xlabel("x-axis").set_ylabel("y-axis").set_zlabel("z-axis");
		g4.plot_equation3d("x*x+y*y");

		Gnuplot g5("lines");
		std::cout << "window 5: splot with hidden3d" << std::endl;
		g5.set_isosamples(25).set_hidden3d();
		g5.plot_equation3d("x*y*y");

		Gnuplot g6("lines");
		std::cout << "window 6: splot with contour" << std::endl;
		g6.set_isosamples(60).set_contour();
		g6.unset_surface().plot_equation3d("sin(x)*sin(y)+4");

		g6.set_surface().replot();

		Gnuplot g7("lines");
		std::cout << "window 7: set_samples" << std::endl;
		g7.set_xrange(-30, 20).set_samples(40);
		g7.plot_equation("besj0(x)*0.12e1").plot_equation("(x**besj0(x))-2.5");

		g7.set_samples(400).replot();

		Gnuplot g8("filledcurves");
		std::cout << "window 8: filledcurves" << std::endl;
		g8.set_legend("outside right top").set_xrange(-5, 5);
		g8.plot_equation("x*x").plot_equation("-x*x+4");

		//
		// Plot an image
		//
		Gnuplot g9;
		std::cout << "window 9: plot_image" << std::endl;
		const int unsigned uiWidth = 255U;
		const int unsigned uiHeight = 255U;
		g9.set_xrange(0, uiWidth).set_yrange(0, uiHeight).set_cbrange(0, 255);
		g9.cmd("set palette gray");
		unsigned char ucPicBuf[uiWidth*uiHeight];
		// generate a greyscale image
		for (unsigned int uiIndex = 0; uiIndex < uiHeight*uiWidth; uiIndex++)
		{
			ucPicBuf[uiIndex] = static_cast<unsigned char>(uiIndex % 255U);
		}
		g9.plot_image(ucPicBuf, uiWidth, uiHeight, "greyscale");

		g9.set_pointsize(0.6).unset_legend().plot_slope(0.8, 20);

		//
		// manual control
		//
		Gnuplot g10;
		std::cout << "window 10: manual control" << std::endl;
		g10.cmd("set samples 400").cmd("plot abs(x)/2"); // either with cmd()
		g10 << "replot sqrt(x)" << "replot sqrt(-x)";    // or with <<

		wait_for_key();
	}
	catch (GnuplotException &ge)
	{
		std::cout << ge.what() << std::endl;
	}
}