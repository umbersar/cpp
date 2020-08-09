// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include<math.h>
#include<float.h>


// TODO: reference additional headers your program requires here
#include<map>
#include<vector>

#include <boost/preprocessor/stringize.hpp>

#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>
using namespace boost::timer;

#include <Eigen/Dense>
using Eigen::Matrix;
using Eigen::MatrixXd;

extern std::map<double(*)(double), std::string> function_map;

void TaylorSeriesUsingNumericalDerivatives(double(*f)(double), double x0, double x);
void GetMachinePrecision();
void MaclaurinSeriesUsingNumericalDerivatives(double(*f)(double), double x);
void sinxMaclaurinSeriesUsingDerivativesFromTigrometry(double x);
void TaylorSeriesUsingSyntheticSubstitutionDerivatives(std::vector<double> polynomialCoefficients, double x0, double x);

double derivativeDividedByFactorialUsingSyntheticSubstitution(std::vector<double> coefficients, double x, int orderOfDerivative);
double derivative(double(*f)(double), double x0, int orderOfDerivative);
double sinDerivativeFromTrignometry(double x0, int orderOfDerivative);

void IterativeRoot(double(*f)(double), double initValue);
void IterativeRootNewtonRaphson(double(*f)(double), double initValue);

long factorial(int i);
double initialGuess(double(*f)(double));
void printMatrix(double **matrix, int rows, int columns);
double** Create2DimensionalMatrix(int matrixRows, int matrixColumns);
void Init2DimensionalMatrix(double **matrix, int matrixRows, int matrixColumns, double *arr, int sizeOfArray);
void Init2DimensionalMatrix(double **matrix, int matrixRows, int matrixColumns, double seed);

double Determinant(double **matrix, int rows, int columns);
void rowReduction(double **matrix, int rowIndex, int rows, int columns);
double** viewOfMatrix(double **matrix, int rows, int columns, int col);
double PivotalCondensationMethod(double **matrix, int rows, int columns);