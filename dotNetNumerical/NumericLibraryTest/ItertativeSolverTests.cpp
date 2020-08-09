#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\NumericLibrary\NumericLibrary.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

//unittesting as explained here:https://msdn.microsoft.com/en-us/library/hh598953.aspx
//Writing Unit tests for C/C++ with the Microsoft Unit Testing Framework for C++

namespace NumericLibraryTest
{
	TEST_CLASS(ItertativeSolverTests)
	{
		//4sinx-e^x=0
		//for Newton Raphson method we DO NOT HAVE TO RESTRUCTURE the equation to x on the LHS.
		//FourSinxMinusexRestructured solved by successive approximation iterative solver should give same result as FourSinxMinusex solved by Newton Raphson iterative solver
		static double FourSinxMinusex(double x) {
			return 4 * sin(x) - exp(x);
		}

		//4sinx-e^x=0
		//solve for x=sin^-1(e^x/4)
		//for iterative solver using successive approximation we have to rewrite the equation to get x on the LHS
		static double FourSinxMinusexRestructured(double x) {
			return asin(exp(x) / 4);
		}

		//x^3+12.1x^2+13.1x+22.2=0
		//for Newton Raphson method we DO NOT HAVE TO RESTRUCTURE the equation to x on the LHS.
		//FourSinxMinusexRestructured solved by successive approximation iterative solver should give same result as FourSinxMinusex solved by Newton Raphson iterative solver
		static double xCubeEquation(double x) {
			return pow(x, 3) + 12.1*pow(x, 2) + 13.1*x + 22.2;
		}

	public:

		TEST_METHOD(IterativeRootSuccessiveApporximationSolverTest)
		{
			//initial value is a guess and generally it is safe to use 0 as initial value
			//this is a iterative solver using successive approximation
			//more efficient iterative solvers also exist like the newton-raphson method
			double initialValue = 0;

			Assert::AreEqual(0.370557, IterativeRoot(FourSinxMinusexRestructured, initialValue), .000001);
		}

		TEST_METHOD(IterativeRootNewtonRaphsonSolverTest)
		{
			//initial value is a guess and generally it is safe to use 0 as initial value
			//for Newton Raphson method we DO NOT HAVE TO RESTRUCTURE the equation to x on the LHS.
			//FourSinxMinusexRestructured solved by successive approximation iterative solver should give same result as FourSinxMinusex solved by Newton Raphson iterative solver
			//note that the difference between Newton Raphson methods and successive approximation is that for successive approximation we restructure the function equation itself to get x on LHS
			//whereas for Newton Raphson method, we modify the taylor series for the function to get x on the LHS

			double initialValue = 0;

			//this is what i got using forward difference approximation in derivative method
			//Assert::AreEqual(0.370557, IterativeRootNewtonRaphson(FourSinxMinusex, initialValue), .000001);
			
			//this is what i got using centered difference approximation(lower order approximation) in derivative method
			Assert::AreEqual(0.370558, IterativeRootNewtonRaphson(FourSinxMinusex, initialValue), .000001);
		}

		TEST_METHOD(IterativeRootSolverComparisonTest)
		{
			double initialValue = 0;

			//Assert::AreEqual(0.370557, IterativeRootNewtonRaphson(FourSinxMinusex, initialValue), .000001);
			//Assert::AreEqual(0.370557, IterativeRoot(FourSinxMinusexRestructured, initialValue), .000001);

			Assert::AreEqual(0.370558, IterativeRootNewtonRaphson(FourSinxMinusex, initialValue), .000001);
			Assert::AreEqual(0.370557, IterativeRoot(FourSinxMinusexRestructured, initialValue), .000001);
		}

		TEST_METHOD(IterativeRootSolverWrongInitialValueComparisonTest)
		{
			double initialValue = 0;

			//0 as the initial guess wont each time. So here for x^3+12.1x^2+13.1x+22.2=0, iterative root could not be found using initial value 0
			//but using -11, we found the root. There are ways to make educated guesses for initial values. page 299 of C.Xavier book on numerical methods explains that.
			Assert::IsTrue(isnan(IterativeRootNewtonRaphson(xCubeEquation, initialValue)));

			initialValue = -11;
			Assert::AreEqual(-11.1, IterativeRootNewtonRaphson(xCubeEquation, initialValue), .000001);
		}

		TEST_METHOD(IterativeRootSolverComputedInitialValueComparisonTest)
		{
			double initialValue = 0;

			//0 as the initial guess wont each time. So here for x^3+12.1x^2+13.1x+22.2=0, iterative root could not be found using initial value 0
			//There are ways to make educated guesses for initial values. page 299 of C.Xavier book on numerical methods explains that.
			initialValue = initialGuess(xCubeEquation);
			Assert::AreEqual(-11.1, IterativeRootNewtonRaphson(xCubeEquation, initialValue), .000001);
		}

	};
}