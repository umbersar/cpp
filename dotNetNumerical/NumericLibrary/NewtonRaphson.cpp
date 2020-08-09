#include "stdafx.h"
#include "NumericLibrary.h"

//given algebraic eq. f(x)=x^7-2x^4+x-8=0. Find it's root. 
//Root of a equation is one which satisfies the equation. To find the root, rewrite the equation so you have x on the LHS. You can rewrite the eq to get x on LHS by many ways.
//If the solution does not converge to the root of the equation, consider rewritting the equation in a different way.
//so you have x=2x^4-x+8-x^7. This is of the form x=f(x0). You use a value of x to find the next value of x by applying the intital guess to the RHS of equation.
//Now successive approximation tells you that you have to keep solving x until the new value of x differs from previous value by 
//whatever tolerance/accuracy the user has specified. 
//Take a guess for the value for x. Let it be 0. And you get x= 8 for the next iteration and then use 8 to find next value of x and so on.

//Newton-Raphson method is a iterative root solver, just like the successive approximation method but it uses taylor series.
//Let f(x)=0 be a given equation. Taylor series for it is:
//taylor series is f(x) = f(x0) + (f'(x0)/1!)(x- x0) + (f''(x0)/2!) (x- x0)^2 + (f'''(0)/3!) (x- x0)^3
//now lets suppose a value x0 (just like we used 0 as the starting point or temporary root) is a educated guess for the root of the equation which then gives us
//the new improved value of root of the equation as x.
//Before we start using the taylor series in iterations, we approximate it by truncating it at the first derivative: f(x) = f(x0) + (f'(x0)/1!)(x- x0). Keep in mind you can 
//think about it in geometrical terms where the derivative is the slope and you are using the slope to find the next f value given a small increase in x.
//Now if using x0 as the initial guess root we get a the new improved root x, then that means f(x)=0.
//0 = f(x0) + (f'(x0)/1!)(x- x0). Rewrite is as:
//x= x0 - f(x0)/f'(x0). Now that becomes an iterative equation.

//this is a iterative solver using successNewton Raphson method
//note that the difference between Newton Raphson methods and successive approximation is that for successive approximation we restructure the function equation itself to get x on LHS
//whereas for Newton Raphson method, we modify the taylor series for the function to get x on the LHSs
NUMERICLIBRARY_API double IterativeRootNewtonRaphson(double(*f)(double), double initValue) {
	int index = 0;
	double newValue;
	double lastValue = initValue;
	while (true) {
		if (index == 500000) {
			printf("Could not solve the equation using Newton Raphson iterative solver in %d iterations. Relax the max iteration count to see if a solution is possible.\n", index);
			printf("The result computed using Newton Raphson iterative solver in %d iterations is %lf\n", index, newValue);
			return NAN;
		}
		newValue = lastValue - f(lastValue) / derivative(f, initValue, 1);
		if (fabs(newValue - lastValue) <= .000001) {
			printf("Root of the equation using Newton Raphson iterative solver found! It is: %lf\n", newValue);
			return newValue;
		}
		else {
			lastValue = newValue;
		}
		index++;
	}
}