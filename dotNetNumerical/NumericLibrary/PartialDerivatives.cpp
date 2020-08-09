#include "stdafx.h"
#include "NumericLibrary.h"


//we are considering a function f(x) and trying to find the change in f(x) for a small change in x. Called the derivative.
//now if the function had more than on arguments, say f(x,y), we can now calculate change in f(x,y) by making a small change in
//x. f'ₓ(x,y) = f(x+delta,y)-f(x,y)/(x+delta-x) or similarly in y. And that is called partial derivation. The rate of change of 
//the function w.r.t one of the variables it depends on. Here also we can use the centred difference formula(and its
//higher order approximations). Now higher order derivatives of f'ₓ(x,y) are f''ₓₓ(x,y) and f'''ₓₓₓ(x,y) and so on.
//How do we define f''ₓᵧ(x,y) means first take the partial derivative w.r.t y and then w.r.t x
//and f''ᵧₓ(x,y) and f''ₓₓ(x,y).



//     |---------|---------|---------|---------|---------|-               
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |---------|---------|---------|---------|---------|-               
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |         | fᵢ₋₁,j₊₁| fᵢ,j₊₁  | fᵢ₊₁,j₊₁|         |        
//     |---------●---------●---------●---------|--------|-               
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |         | fᵢ₋₁,j  | fᵢ,j    | fᵢ₊₁,j  |         |        
//     |---------●---------●---------●---------|--------|-               
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |         | fᵢ₋₁,j₋₁|fᵢ,j₋₁   | fᵢ₊₁,j₋₁|         |        
//     |---------●---------●---------●---------|--------|-               
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |         |         |         |         |         |        
//     |---------|---------|---------|---------|---------|-               


//Finite difference approximations of partial derivatives using central differences:
//Conside a mesh of points in xy plane with Δx=h and Δy=k. 
//Uₓ = f'ₓ(x,y) = (f(x+h,y)-f(x-h,y))/(2*h) = (fᵢ₊₁,j - fᵢ₋₁,j) / 2*h 
//Uₓₓ = f''ₓₓ(x,y) = (f(x+h,y) -2*f(x,y) + f(x-h,y))/(2*h) = (fᵢ₊₁,j - 2*fᵢ,j + fᵢ₋₁,j) / 2*h
//
//Uᵧ = f'ᵧ(x,y) = (f(x,y+k)-f(x,y-k))/(2*k) =  (fᵢ,j₊₁ - fᵢ,j₋₁) / 2*k
//Uᵧᵧ = f''ᵧᵧ(x,y) = (f(x+h,y) -2*f(x,y) + f(x-h,y))/(2*h) = (fᵢ,j₊₁ - 2*fᵢ,j + fᵢ,j₋₁) / 2*k
//
//Uₓᵧ = f''ₓᵧ(x,y) = (f(x+h,y+k) - f(x+h,y-k) - f(x-h,y+k) + f(x-h,y-k))/(4*h*k) = (fᵢ₊₁,j₊₁ - fᵢ₊₁,j₋₁ - fᵢ₋₁,j₊₁ + fᵢ₋₁,j₋₁) / (4*h*k) 
//
//The value Uₓₓ + Uᵧᵧ is used for a number of problems and is called laplacian of U. It is denoted by ∇².
//∇²U = Uₓₓ + Uᵧᵧ 
//laplace equation is ∇²U = Uₓₓ + Uᵧᵧ = 0. Note the difference between laplace operator(laplacian) and laplace equation. 
//Todo: there is a connection of laplace operator with eigen values.

//let f(x,y)= x^2 +y 
//f'ₓ(x,y) = f(x+delta,y)-f(x,y)/(x+delta-x).
//1. For f'ₓ(4,5) the call is partialDerivative(f, {4,5}, {0}, 1 , 1)
//2. For f'ₓᵧ(4,5) the call is partialDerivative(f, {4,5}, {1,0}, 2 , 1).
// I assume that the derivative w.r.t y is calculated before w.r.t x. Otherwise the call would have been partialDerivative(f, {4,5}, {0,1}, 2 , 1)
// Assuming the first is true for order in which the derivatives are computed, then how is the return value of first derivative used? 
// Because after we compute the first partial derivative, the function f gets evaluated to value and we loose the function.
// Finite differences give you partial derivative approximations formulas as mentioned here https://en.wikipedia.org/wiki/Finite_difference#Finite_difference_in_several_variables
//3. For f''ₓᵧ(4,5) the call is partialDerivative(f, {4,5}, {1,0}, 2 , 2)


NUMERICLIBRARY_API double partialDerivative(double(*f)(double, double), double *varValues, /*int sizeOfVarValues,*/ int *changeIndex,
	int sizeOfchangeIndex, int orderOfDerivative) {
	throw NotImplementedException("partial derivatives not implmented yet.");
}
//
//NUMERICLIBRARY_API double partialDerivative(double(*f)(double, double), double *varValues, /*int sizeOfVarValues,*/ int *changeIndex,
//	int sizeOfchangeIndex, int orderOfDerivative) {
//	//how to choose meaningful delta 
//	double delta = 0.0625;
//	int sizeOfvarValues = 2;
//
//	if (orderOfDerivative == 1) {// && sizeOfchangeIndex == 1) {
//		//change the var to be w.r.t which the derivative is being computed. And pass in the rest of the vars as such
//		//you have to know the position of the var being hange in the argument list as well. 
//		double *terms = (double *)malloc(sizeOfvarValues * sizeof(double));
//		int varToChangeIndex = -1;
//		for (int i = 0; i < sizeOfchangeIndex; i++) {
//			for (int j = 0; j < sizeOfvarValues; j++) {
//				if (j == changeIndex[i]) {
//					terms[j] = varValues[j] + delta;
//					varToChangeIndex = j;
//				}
//				else
//					terms[j] = varValues[j];
//			}
//		}
//		return (terms[0] - terms[1]) / ((varValues[varToChangeIndex] + delta) - varValues[varToChangeIndex]);
//	}
//	else {
//		return (derivative(f, x0 + delta, orderOfDerivative - 1) - derivative(f, x0, orderOfDerivative - 1))
//			/ (x0 + delta - x0);
//	}
//}

