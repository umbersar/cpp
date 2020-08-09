#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\NumericLibrary\NumericLibrary.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

//unittesting as explained here:https://msdn.microsoft.com/en-us/library/hh598953.aspx
//Writing Unit tests for C/C++ with the Microsoft Unit Testing Framework for C++

namespace NumericLibraryTest
{
	TEST_CLASS(DifferentialEquationsTests)
	{
		static double xPlusY(double x, double y)
		{
			return x + y;
		}

		static double xPlusYSqr(double x, double y)
		{
			return x + y*y;
		}

		static double firstDifferentialEq(double x, double y, double z)
		{
			return z;
		}

		static double secondDifferentialEq(double x, double y, double z)
		{
			return 6 * y - z;
		}

		static double u1DifferentialEq(double t, double u1, double u2)
		{
			return -4 * u1 - 2 * u2 + cos(t) + 4 * sin(t);
		}

		static double u2DifferentialEq(double t, double u1, double u2)
		{
			return 3 * u1 + u2 - 3 * sin(t);
		}

		static	double α;
		static	double β;
		static	double δ;
		static	double Ɣ;

		//α, β, δ and Ɣ constants have been scoped to simulate closures
		static double rabbitsDifferentialEq(double t, double x, double y)//, double α, double β)
		{
			return α*x - β*x*y;
		}

		//α, β, δ and Ɣ constants have been scoped to simulate closures
		static double foxesDifferentialEq(double t, double x, double y)//, double δ, double Ɣ)
		{
			return δ*x*y - Ɣ*y;
		}

	public:

		TEST_METHOD(EulersInitialValueMethodTest)
		{
			//the differential equation is dy/dx = x+y
			//initial conditions are 
			double x0 = 0;
			double y0 = 1;

			//we have to find value for
			double x = 1;
			double y = EulersMethod(x0, y0, x, xPlusY);

			//the exact solution is 3.44
			Assert::AreEqual(3.43656, y, .00001);
			Assert::AreNotEqual(3.44, y);

		}

		TEST_METHOD(EulersInitialValueModifiedMethodTest)
		{
			//the differential equation is dy/dx = x+y
			//initial conditions are 
			double x0 = 0;
			double y0 = 1;

			//we have to find value for
			double x = 1;
			double y = EulersMethodModified(x0, y0, x, xPlusY);

			//the exact solution is 3.44
			Assert::AreEqual(3.44107, y, .00001);
			Assert::AreEqual(3.44, y, .01);

		}

		TEST_METHOD(RungeKuttaTest)
		{
			//the differential equation is dy/dx = x+y^2
			//initial conditions are 
			double x0 = 0;
			double y0 = 1;

			//we have to find value for
			double x = .1;
			double y = RungeKuttaFourthOrder(x0, y0, x, xPlusYSqr);

			Assert::AreEqual(1.1165, y, .0001);

			x0 = .1;
			y0 = y;
			x = .2;
			y = RungeKuttaFourthOrder(x0, y0, x, xPlusYSqr);
			Assert::AreEqual(1.2736, y, .0001);

			x0 = 0;
			y0 = 1;
			x = .2;
			y = RungeKuttaFourthOrder(x0, y0, x, xPlusYSqr);
			Assert::AreEqual(1.2736, y, .0001);
		}

		TEST_METHOD(RungeKuttaForTwoDifferentialEquationsTest)
		{
			//we have two first order differential equations: dy/dx=z and dz/dx=6y−z. Inititial conditions are given
			//these 2 first order differential equations were initially a single second order equation: d²y/dx² + dy/dx −6y=0
			//so runge kutta requires the equations to be in first order.

			//https://www.math.ohiou.edu/courses/math3600/lecture29.pdf https://math.berkeley.edu/~zworski/128/psol12.pdf (in git at rungekutta1.pdf and rungekutta2.pdf)

			//initial conditions are
			double x0 = 0;
			double y0 = 3;
			double z0 = 1;

			//we have to find value for
			double x = 1;

			//define h as well as malloc the variable arrays accordingly
			double h = .1;
			int sizeOfVarArrays = ceil((x - x0) / h) + 1;
			double *xVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			double *yVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			double *zVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			xVarArray[0] = x0;
			yVarArray[0] = y0;
			zVarArray[0] = z0;
			xVarArray[sizeOfVarArrays - 1] = x;

			double(*f[2])(double, double, double) = { firstDifferentialEq,secondDifferentialEq };
			RungeKuttaForTwoDifferentialEquationsFourthOrder(xVarArray, yVarArray, zVarArray, sizeOfVarArrays, f, h);

			Assert::AreEqual(14.82757, yVarArray[sizeOfVarArrays - 1], .00001);

			free(xVarArray);
			free(yVarArray);
			free(zVarArray);
		}

		TEST_METHOD(RungeKuttaForTwoDifferentialEquationsWorkedExampleTest)
		{
			//worked example is avaiable here https://math.berkeley.edu/~zworski/128/psol12.pdf

			//initial conditions are
			double t0 = 0;
			double u10 = 0;
			double u20 = -1;

			//we have to find value for
			double t = 2;

			//define h as well as malloc the variable arrays accordingly
			double h = .1;
			int sizeOfVarArrays = ceil((t - t0) / h) + 1;
			double *tVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			double *u1VarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			double *u2VarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			tVarArray[0] = t0;
			u1VarArray[0] = u10;
			u2VarArray[0] = u20;
			tVarArray[sizeOfVarArrays - 1] = t;

			double(*f[2])(double, double, double) = { u1DifferentialEq,u2DifferentialEq };
			RungeKuttaForTwoDifferentialEquationsFourthOrder(tVarArray, u1VarArray, u2VarArray, sizeOfVarArrays, f, h);

			Assert::AreEqual(1.14332436, u1VarArray[sizeOfVarArrays - 1], .00000001);
			Assert::AreEqual(-0.36936318, u2VarArray[sizeOfVarArrays - 1], .00000008);

			free(tVarArray);
			free(u1VarArray);
			free(u2VarArray);
		}

		TEST_METHOD(RungeKuttaForPredatorPreyEquationsTest)
		{
			//x is number of rabbits and y is number of foxes and t represents time.
			//equation for rabbits is: dx/dt = αx -βxy gives the rate of change of population of x(rabbits)
			//equation for foxes is: dy/dt = δxy -Ɣy gives the rate of change of population of y(foxes)
			//where α, β, δ and Ɣ are constants describing the interation of the 2 species

			//pick some values for the constants...todo: try out the closures to capture these constants and then you
			//would get the 2 functions which only need the 2 argumest x,y and t to be passed and not the constants
			//https://www.youtube.com/watch?v=0LzDiScAcJI&index=12&list=PLpqe4qvXoGHQaClB8ImfPno8vN5R89RIQ
			// instead of closure, lets use global variable scoping to give similar effect. http://programmers.stackexchange.com/questions/263057/every-function-is-a-closure

			/*double α = 2;
			double β = 2;
			double δ = 1;
			double Ɣ = 1;*/

			//initial conditions are 3 rabbits, 2 foxes at time 0
			double t0 = 0;
			double x0 = 3;
			double y0 = 2;

			//we have to find value for rabbits and foxes at t=.1
			double t = .1;
			//define h as well as malloc the variable arrays accordingly
			double h = .1;
			int sizeOfVarArrays = ceil((t - t0) / h) + 1;
			double *tVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			double *xVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			double *yVarArray = (double*)malloc(sizeOfVarArrays * sizeof(double));
			tVarArray[0] = t0;
			xVarArray[0] = x0;
			yVarArray[0] = y0;
			tVarArray[sizeOfVarArrays - 1] = t;
			double(*f[2])(double, double, double) = { rabbitsDifferentialEq,foxesDifferentialEq };
			RungeKuttaForTwoDifferentialEquationsFourthOrder(tVarArray, xVarArray, yVarArray, sizeOfVarArrays, f, h);

			Assert::AreEqual(2.364, xVarArray[sizeOfVarArrays - 1], .001);
			Assert::AreEqual(2.367, yVarArray[sizeOfVarArrays - 1], .001);


			free(xVarArray);
			free(yVarArray);
			free(tVarArray);
		}

		TEST_METHOD(AdamsBashforthMethodTest)
		{
			//http://math.stackexchange.com/a/721234/377760
			//a single second order equation : d²y / dx² + dy / dx −6y = 0 gets decomposed into two first order differential 
			//equations: dy/dx=z and dz/dx=6y−z. Because runge kutta requires the equations to be in first order.

			//https://www.math.ohiou.edu/courses/math3600/lecture29.pdf https://math.berkeley.edu/~zworski/128/psol12.pdf (in git at rungekutta1.pdf and rungekutta2.pdf)

			//initial conditions are
			double x0 = 0;
			double y0 = 3;
			double z0 = 1;

			//we have to find value for
			double x = 1;

			//define h as well as malloc the variable arrays accordingly
			double h = .1;
			double *xVarArray;// = 0;
			double *yVarArray;// = 0;
			double *zVarArray;// = 0;
			int *sizeOfVarArrays;// = 0;s

			double(*f[2])(double, double, double) = { firstDifferentialEq,secondDifferentialEq };
			AdamsBashforthMoultonFourthOrderForTwoDifferentialEquations(x0, y0, z0, x, &xVarArray, &yVarArray, &zVarArray, &sizeOfVarArrays, f, h);

			Assert::AreEqual(14.82889, yVarArray[*sizeOfVarArrays - 1], .00001);

			free(sizeOfVarArrays);
			free(xVarArray);
			free(yVarArray);
			free(zVarArray);
		}
	};
	double DifferentialEquationsTests::α = 2;
	double DifferentialEquationsTests::β = 2;
	double DifferentialEquationsTests::δ = 1;
	double DifferentialEquationsTests::Ɣ = 1;
}