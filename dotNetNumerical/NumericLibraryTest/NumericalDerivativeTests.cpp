#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\NumericLibrary\NumericLibrary.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

//unittesting as explained here:https://msdn.microsoft.com/en-us/library/hh598953.aspx
//Writing Unit tests for C/C++ with the Microsoft Unit Testing Framework for C++

namespace NumericLibraryTest
{
	TEST_CLASS(NumericalDerivativeTests)
	{
		static double negativeExp(double x)
		{
			return exp(-x);
		}

		//y=(1-x^2)^.5
		static double OneMinusXSq(double x)
		{
			return pow((1 - pow(x, 2)), .5);
		}

		//y=1/(1+x)^.5
		static double OnePlusXSqrt(double x)
		{
			return 1 / pow(1 + x, .5);
		}

		//polynomial f(x) = 3x ^ 3 + 2x ^ 2 + x - 2
		static double polynomialFunc(double x)
		{
			return 3 * pow(x, 3) + 2 * pow(x, 2) + x - 2;
		}

	public:

		TEST_METHOD(MaclaurinSeriesUsingNumericalDerivativesTest)
		{
			//this is what i got using forward difference approximation. Using derivative method
			/*Assert::AreEqual(0.245322, MaclaurinSeriesUsingNumericalDerivatives(sin, .25), .000001);
			Assert::AreEqual(0.0684, MaclaurinSeriesUsingNumericalDerivatives(cos, 1.5, .00001, 15), .001);*/

			//this is what i got using centered difference approximation(lower order approximation). Using derivative1 method
			//Assert::AreEqual(0.247404, MaclaurinSeriesUsingNumericalDerivatives(sin, .25), .000001);
			//Assert::AreEqual(0.0717111, MaclaurinSeriesUsingNumericalDerivatives(cos, 1.5, .00001, 15), .000001);

			//this is what i got using centered difference approximation(higher order approximation). Using derivative2 method
			Assert::AreEqual(0.247404, MaclaurinSeriesUsingNumericalDerivatives(sin, .25), .000001);
			Assert::AreEqual(0.070738, MaclaurinSeriesUsingNumericalDerivatives(cos, 1.5, .00001, 15), .000001);
		}

		TEST_METHOD(MaclaurinSeriesUsingDerivativesFromTigrometryTest)
		{
			Assert::AreEqual(0.247404, sinxMaclaurinSeriesUsingDerivativesFromTigrometry(.25), .0000001);
		}

		TEST_METHOD(MaclaurinSeriesUsingNumericalDerivativesForCustomFunctionTest)
		{
			//this is what i got using forward difference approximation. Using derivative method
			/*Assert::AreEqual(0.823758, MaclaurinSeriesUsingNumericalDerivatives(negativeExp, .2), .000001);
			Assert::AreEqual(0.836558, MaclaurinSeriesUsingNumericalDerivatives(OneMinusXSq, .5, .0001, 20), .000001);
			Assert::AreEqual(0.74118, MaclaurinSeriesUsingNumericalDerivatives(OnePlusXSqrt, 1 / pow(1.4, .5), .01, 15), .00001);*/

			//this is what i got using centered difference approximation(lower order approximation). Using derivative1 method
			//Assert::AreEqual(0.818624, MaclaurinSeriesUsingNumericalDerivatives(negativeExp, .2), .000001);
			//Assert::AreEqual(0.864975, MaclaurinSeriesUsingNumericalDerivatives(OneMinusXSq, .5, .0001, 20), .000001);
			//Assert::AreEqual(0.845, MaclaurinSeriesUsingNumericalDerivatives(OnePlusXSqrt, .4, .01, 15), .001);

			//this is what i got using centered difference approximation(higher lower order approximation). Using derivative2 method
			Assert::AreEqual(0.818731, MaclaurinSeriesUsingNumericalDerivatives(negativeExp, .2), .000001);
			Assert::AreEqual(0.8662, MaclaurinSeriesUsingNumericalDerivatives(OneMinusXSq, .5, .001, 11), .0001);
			Assert::AreEqual(0.845, MaclaurinSeriesUsingNumericalDerivatives(OnePlusXSqrt, .4, .01, 15), .001);
		}

		TEST_METHOD(TaylorSeriesUsingNumericalDerivativesForLogFunctionTest)
		{
			//taylor series where x0=1, f(x) is log x and  we have to solve for f(x) where x is 2


			//this is what i got using forward difference approximation. Using derivative method
			//Assert::AreEqual(0.69, TaylorSeriesUsingNumericalDerivatives(log, 1, 2, .01, 15), .01);

			//this is what i got using centered difference approximation(lower order approximation). Using derivative1 method
			//Assert::AreEqual(0.69, TaylorSeriesUsingNumericalDerivatives(log, 1, 2, .01, 15), .01);

			//this is what i got using centered difference approximation(higher lower order approximation). Using derivative2 method
			//cant solve log using centered difference approximation even for less accuracy. forward difference approximation would have worked.
			//and if we had tried to solve it using more more iterations, we would have gotten a divide by 0 error as the x0 in the 
			//derivative calculation method would have gotten smaller and smaller. Check the method ExceptionsInNumericalDerivativesTest
			//down below.
			double val = TaylorSeriesUsingNumericalDerivatives(log, 1, 2, .1, 8);
			Assert::IsTrue(isnan(val));
		}

		TEST_METHOD(ExceptionsInNumericalDerivativesTest)
		{
			//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
			
			wchar_t message[200];
			try
			{
				//this is supposed to fail for the derivative calculation centered difference approximation(higher or lower order approximation)
				//because while computing higher order derivatives, we will reach a point where we get divide by 0 floating error.
				//As x0 will tend to 0, the function log, of which the derivative is being computed, will do a divide by 0.
				//this wont happen were we using forward difference approximation
				TaylorSeriesUsingNumericalDerivatives(log, 1, 2, .01, 12);
				
				_swprintf(message, L"Expected NotComputableException");
				Assert::Fail(message, LINE_INFO());
			}
			catch (NotComputableException ex)
			{
				//swallow the exception
				//ex.what();
			}
			catch (...)
			{
				Assert::Fail(message, LINE_INFO());
			}

			try
			{
				MaclaurinSeriesUsingNumericalDerivatives(OneMinusXSq, .5, .0001, 15);

				_swprintf(message, L"Expected NotComputableException");
				Assert::Fail(message, LINE_INFO());
			}
			catch (NotComputableException ex)
			{
				//swallow the exception
				//ex.what();
			}
			catch (...)
			{
				_swprintf(message, L"Incorrect exception. Expected NotComputableException");
				Assert::Fail(message, LINE_INFO());
			}
		}

		TEST_METHOD(SolvePolynomialUsingSyntheticSubstitutionTest)
		{
			//synthetic substitution can be used for solving polynomials. See the 'computational science' word doc. here we are passing the orderOfDerivative parameter set to 0. 
			//that means just calculate the polynomial.
			std::vector<double> polynomialCoefficients = { 3,2,1,-2 };

			double dsynsub0 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 0);
			printf("The solution of polynomial calculated using synthetic substitution is %lf\n", dsynsub0);
			Assert::AreEqual(32, dsynsub0, .01);

			dsynsub0 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 3, 0);
			printf("The solution of polynomial calculated using synthetic substitution is %lf\n", dsynsub0);
			Assert::AreEqual(100, dsynsub0, .01);
		}

		TEST_METHOD(ComputeDerivativeUsingSyntheticSubstitutionTest)
		{
			//if the function is a polynomial, then instead of numerical derivaltives, we can also use Synthetic Substitution Derivatives
			//this catually gives us the derivative/factorial. So to get the exact derivative value, you need to, multiply with factorial value
			std::vector<double> polynomialCoefficients = { 3,2,1,-2 };

			double dsynsub1 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 1) * factorial(1);
			Assert::AreEqual(45, dsynsub1, .01);

			double dsynsub2 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 2) * factorial(2);
			Assert::AreEqual(40, dsynsub2, .01);

			double dsynsub3 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 3) * factorial(3);
			Assert::AreEqual(18, dsynsub3, .01);
		}

		TEST_METHOD(CompareDerivativeComputedUsingSyntheticSubstitutionAndNumericalMethodsTest)
		{
			//if the function is a polynomial, then instead of numerical derivaltives, we can also use Synthetic Substitution Derivatives
			//this catually gives us the derivative/factorial. So to get the exact derivative value, you need to, multiply with factorial value
			std::vector<double> polynomialCoefficients = { 3,2,1,-2 };

			double dsynsub1 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 1) * factorial(1);
			Assert::AreEqual(45, dsynsub1, .01);

			double dsynsub2 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 2) * factorial(2);
			Assert::AreEqual(40, dsynsub2, .01);

			double dsynsub3 = derivativeDividedByFactorialUsingSyntheticSubstitution(polynomialCoefficients, 2, 3) * factorial(3);
			Assert::AreEqual(18, dsynsub3, .01);

			double dnum1 = derivative(polynomialFunc, 2, 1);
			double dnum2 = derivative(polynomialFunc, 2, 2);
			double dnum3 = derivative(polynomialFunc, 2, 3);

			//this is what i got using forward difference approximation. Using derivative method
			//Assert::AreEqual(46.261719, dnum1, .000001);
			//Assert::AreEqual(41.125000, dnum2, .000001);
			//Assert::AreEqual(18.00, dnum3, .000001);

			//this is what i got using centered difference approximation(lower order approximation). Using derivative1 method
			/*Assert::AreEqual(45.0117, dnum1, .0001);
			Assert::AreEqual(40, dnum2, .000001);
			Assert::AreEqual(18.00, dnum3, .000001);*/

			//this is what i got using centered difference approximation(higher lower order approximation). Using derivative2 method
			Assert::AreEqual(45, dnum1, .0001);
			Assert::AreEqual(40, dnum2, .000001);
			Assert::AreEqual(18.00, dnum3, .000001);

			printf("The derivatives calculated using synthetic substitution are %lf, %lf and %lf\n", dsynsub1, dsynsub2, dsynsub3);
			printf("The derivatives calculated using numerical methods are %lf, %lf and %lf\n", dnum1, dnum2, dnum3);


		}

		TEST_METHOD(CompareTaylorSeriesUsingSyntheticSubstitutionDerivativesAndNumericalDerivativesTest)
		{
			//taylor series where x0=2, f(x) is a polynomial f(x) = 3x^3 + 2x^2 + x - 2 
			//and we have to solve for f(x) where x is 3
			std::vector<double> polynomialCoefficients = { 3,2,1,-2 };

			double taySSS = TaylorSeriesUsingSyntheticSubstitutionDerivatives(polynomialCoefficients, 2, 3);
			Assert::AreEqual(100, taySSS, .000001);


			//taylor series where x0=2, f(x) is a polynomial f(x) = 3x^3 + 2x^2 + x - 2 and  we have to solve for f(x) where x is 3
			//polynomialFunc and polynomialCoefficients denote same polynomial
			double taySNM = TaylorSeriesUsingNumericalDerivatives(polynomialFunc, 2, 3);
			//this is what i got using forward difference approximation. Using derivative method
			//Assert::AreEqual(101.824219, taySNM, .000001);
			//this is what i got using centered difference approximation(lower order approximation). Using derivative1 method
			//Assert::AreEqual(100.012, taySNM, .001);
			//this is what i got using centered difference approximation(higher order approximation). Using derivative2 method
			Assert::AreEqual(100, taySNM, .001);

			//similarly calculate MaclaurinSeriesUsingNumericalDerivatives and MaclaurinSeriesUsingSyntheticSubstitution

		}

	};
}