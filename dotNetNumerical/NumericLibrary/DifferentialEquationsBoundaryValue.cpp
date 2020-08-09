#include "stdafx.h"
#include "NumericLibrary.h"

//differential equations: they can be ODE(first order or higher order) or partial(first order or higher order). 

//Initial Value problem: The differential equation has to be first order differential equation. If it is not, then convert
//the higher order differential equation to a system of siumlataneous first order differential equations. And then solve
//the system of simulatenous equations using any of the initial value diff. eq. solvers
//fᵢ₊₁(x) = fᵢ(x) + f'ᵢ(x)*h. Given initial values of both x and y and the differential equation f'ᵢ(x), we can select a small h
//and keep iterating over the equation to plot the function as well as find the value of y at a particualar point
//So how do we solve intial value partial differential equations(PDEs)?

//Boundary value problem(BVP): fᵢ₊₁(x) = fᵢ(x) + f'ᵢ(x)*h. In this case the differential equation, an inital conditions x0 and y0 and
//a terminal conditions xn and yn are given and we can select a small h to divide the range xn to x0 in small intervals. Now in 
//this case the differential equation can either be ODE or PDE. Here we consider an ODE and specifically 2nd order ODE(what about 
//first order or higher order ODE?).
//
//One popular way to solve 2nd order ODE BVP is to use finite difference method. What we do is that we replace the derivatives occuring 
//in the differential equation and the boundary conditions is replace by their finite difference approximations and the resulting system
//of linear equations is solved by any standard procedure.
//refer NumericalFiniteDifferenceDerivatives.cpp
//Now the formulas for central difference approximations to derivatives are:
//f'(x0)=(f(x0 + delta) - f(x0 - delta))/ 2delta 
//f''(x0)=(f(x0 + delta) - 2*f(x0) + f(x0 - delta)) / (delta^2)
//f'''(x0)=(f(x0 + 2*delta) - 2*f(x0 + delta) + 2*f(x0 - delta) - f(x0 - 2*delta)) / (2*delta^3)
//f'''(x0)=(f(x0 + 2*delta) - 4*f(x0 + delta) + 6*f(x0) - 4*f(x0 - *delta) + f(x0 - 2*delta)) / (delta^4)
//writting these central differences approximations to derivatives in terms of index we get:
//f'ᵢ(x)=(fᵢ₊₁(x) - fᵢ₋₁(x)) / 2*delta 
//f''ᵢ(x)=(fᵢ₊₁(x) - 2 * fᵢ(x) + fᵢ₋₁(x)) / delta ^ 2 
//f'''ᵢ(x)=(fᵢ₊₂(x) - 2 * fᵢ₊₁(x) + 2*fᵢ₋₁(x) - fᵢ₋₂(x)) / 2*delta^3 
//f''''ᵢ(x)=(fᵢ₊₂(x) - 4 * fᵢ₊₁(x) + 6*fᵢ(x) - 4*fᵢ₋₁(x) + fᵢ₋₂(x)) / delta^4
//
//Example: solve y'' = x + y with bounday conditions y(0)=y(1)=0
//Lets assume h=.25. Use the central difference approximation in place of the derivative in the differential eq
//(fᵢ₊₁(x) - 2*fᵢ(x) + fᵢ₋₁(x)) / delta^2  = xᵢ + fᵢ(x)
//(yᵢ₊₁ - 2*yᵢ + yᵢ₋₁) / delta^2  = xᵢ + yᵢ
//(yᵢ₊₁ - 2*yᵢ + yᵢ₋₁) / delta^2  - yᵢ= xᵢ 
//now use the step size value, it becomes
//16yᵢ₊₁ - 33*yᵢ + 16yᵢ₋₁ = xᵢ 
//
//Due to the step size we have choosen, the variables would be:
//x0,y0, x1,y1, x2,y2, x3,y3 and x4,y4 and we would have 4 equations as shown below :
//16yᵢ₊₁ - 33*yᵢ + 16yᵢ₋₁ = xᵢ where i = 0 to 4
//16y₂ - 33*0 + 16y₋₁ = 0 for i=0. Now since 16y₋₁ refers to a solution outside the boundary, we would disregard this equation.
//
//16y₂ - 33*y₁ + 16(0) = .25 for i=1.
//16y₃ - 33*y₂ + 16y₁ = .50 for i=2.
//16(0) - 33*y₃ + 16y₂ = .75 for i = 3.
//
//16y₅ - 33*y₄ + 16y₃ = 1 for i = 4. Now since 16y₅ refers to a solution outside the boundary, we would disregard this equation.
//So we are left with 3 linear equations for the system(for i=1 to 3):
//we can solve it, for example, using LUDecompositionLinearEquationsSolver coded in MatrixLinearAlgebra.cpp
//Sample usage for the solver can be seen in test LULinearEquationSolverTest
