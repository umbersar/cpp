#include "stdafx.h"
#include "NumericLibrary.h"

//given algebraic eq. f(x)=x^7-2x^4+x-8=0. Find it's root. 
//Root of a equation is one which satisfies the equation. To find the root, rewrite the equation so you have x on the LHS. You can rewrite the eq to get x on LHS by many ways.
//If the solution does not converge to the root of the equation, consider rewritting the equation in a different way.
//so you have x=2x^4-x+8-x^7. This is of the form x=f(x0). You use a value of x to find the next value of x by applying the intital guess to the RHS of equation.
//Now successive approximation tells you that you have to keep solving x until the new value of x differs from previous value by 
//whatever tolerance/accuracy the user has specified. 
//Take a guess for the value for x. Let it be 0. And you get x= 8 for the next iteration and then use 8 to find next value of x and so on.


//this is a iterative solver using successive approximation
//more efficient iterative solvers also exist like the newton-raphson method, Chord's method and bisection method
//if f(x) evaluation is easy, then use newton raphson. Otherwise  use chord's method and verify convergence using bisection method once every 10 or so iteration.
//bisection method is the most reliable one in that it will certainly find the root but is very slow.
NUMERICLIBRARY_API double IterativeRoot(double(*f)(double), double initValue) {
	int index = 0;
	double newValue;
	double lastValue = initValue;
	while (true) {
		if (index == 500000) {
			printf("Could not solve the equation using iterative solver in %d iterations. Relax the max iteration count to see if a solution is possible\n", index);
			return NAN;
		}
		newValue = f(lastValue);
		if (fabs(newValue - lastValue) <= .000001) {
			printf("Root of the equation using iterative solver found! It is: %lf\n", newValue);
			return newValue;
		}
		else {
			lastValue = newValue;
		}
		index++;
	}
}

