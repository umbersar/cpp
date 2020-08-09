// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the NUMERICLIBRARY_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// NUMERICLIBRARY_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef NUMERICLIBRARY_EXPORTS
#define NUMERICLIBRARY_API __declspec(dllexport)
#else
#define NUMERICLIBRARY_API __declspec(dllimport)
#endif

// This class is exported from the NumericLibrary.dll
class NUMERICLIBRARY_API CNumericLibrary {
public:
	CNumericLibrary(void);
	// TODO: add your methods here.
};

extern NUMERICLIBRARY_API int nNumericLibrary;

NUMERICLIBRARY_API int fnNumericLibrary(void);


class NUMERICLIBRARY_API NotImplementedException :public std::logic_error {
public:
	NotImplementedException(const char * message);
};

class NUMERICLIBRARY_API NotComputableException :public std::runtime_error {
public:
	NotComputableException(const char * message);
};

NUMERICLIBRARY_API double TaylorSeriesUsingNumericalDerivatives(double(*f)(double), double x0, double x, double accuracy = .00001, int maxIterations = 10);
NUMERICLIBRARY_API void GetMachinePrecision();
NUMERICLIBRARY_API double MaclaurinSeriesUsingNumericalDerivatives(double(*f)(double), double x, double accuracy = .00001, int maxIterations = 10);
NUMERICLIBRARY_API double sinxMaclaurinSeriesUsingDerivativesFromTigrometry(double x, double accuracy = .00001, int maxIterations = 10);
NUMERICLIBRARY_API double TaylorSeriesUsingSyntheticSubstitutionDerivatives(std::vector<double> polynomialCoefficients, double x0, double x, double accuracy = .00001, int maxIterations = 10);

NUMERICLIBRARY_API double derivativeDividedByFactorialUsingSyntheticSubstitution(std::vector<double> coefficients, double x, int orderOfDerivative);
NUMERICLIBRARY_API double derivative(double(*f)(double), double x0, int orderOfDerivative);
NUMERICLIBRARY_API double sinDerivativeFromTrignometry(double x0, int orderOfDerivative);

NUMERICLIBRARY_API double IterativeRoot(double(*f)(double), double initValue);
NUMERICLIBRARY_API double IterativeRootNewtonRaphson(double(*f)(double), double initValue);

NUMERICLIBRARY_API long factorial(int i);
NUMERICLIBRARY_API double initialGuess(double(*f)(double));
NUMERICLIBRARY_API void printMatrix(double **matrix, int rows, int columns);
NUMERICLIBRARY_API double** Create2DimensionalMatrix(int matrixRows, int matrixColumns);
NUMERICLIBRARY_API void Init2DimensionalMatrix(double **matrix, int matrixRows, int matrixColumns, double *arr, int sizeOfArray);
NUMERICLIBRARY_API void Init2DimensionalMatrix(double **matrix, int matrixRows, int matrixColumns, double seed);

NUMERICLIBRARY_API enum TriangularForm { Upper, Lower, Diagonal };
NUMERICLIBRARY_API double Determinant(double **matrix, int rows, int columns);
NUMERICLIBRARY_API void rowReduction(double **matrix, int rowIndex, int rows, int columns, TriangularForm tf, bool usePivoting, double *rhs = nullptr);
NUMERICLIBRARY_API void reduceToTriangularForm(double **matrix, int rows, int columns, TriangularForm tf, bool usePivoting, double *rhs = nullptr);
NUMERICLIBRARY_API void reduceToIdentityMatrix(double **matrix, int rows, int columns, double *variables, double *rhs);
NUMERICLIBRARY_API void GaussElimination(double **matrix, int rows, int columns, double *variables, bool usePivoting, double *rhs);
NUMERICLIBRARY_API void GaussJordanElimination(double **matrix, int rows, int columns, double *variables, bool usePivoting, double *rhs);
NUMERICLIBRARY_API double** viewOfMatrix(double **matrix, int rows, int columns, int col);
NUMERICLIBRARY_API double PivotalCondensationMethod(double **matrix, int rows, int columns);
NUMERICLIBRARY_API bool isLowerTriangularMatrix(double **matrix, int rows, int columns);
NUMERICLIBRARY_API bool isUpperTriangularMatrix(double **matrix, int rows, int columns);
NUMERICLIBRARY_API bool isDiagonalMatrix(double **matrix, int rows, int columns);
NUMERICLIBRARY_API void ForwardSubstitution(double **matrix, int rows, int columns, double *variables, double *rhs);
NUMERICLIBRARY_API void BackwardSubstitution(double **matrix, int rows, int columns, double *variables, double *rhs);
NUMERICLIBRARY_API double** InverseOfMatrixUsingGaussElimination(double **matrix, int rows, int columns);
NUMERICLIBRARY_API void rowReduction(double **matrix, int rowIndex, int rows, int columns, TriangularForm tf, bool usePivoting, double **rhs);
NUMERICLIBRARY_API void reduceToTriangularForm(double **matrix, int rows, int columns, TriangularForm tf, bool usePivoting, double **rhs);
NUMERICLIBRARY_API void BackwardSubstitution(double **matrix, int rows, int columns, double **inverseOfMatrix, double **rhsIdentityMatrix);
NUMERICLIBRARY_API void LUFactorize(double **matrix, int rows, int columns, double **L, double **U);
NUMERICLIBRARY_API void LUDecompositionLinearEquationsSolver(double **matrix, int rows, int columns, double *variables, double *rhs);
NUMERICLIBRARY_API double LinearInterpolation(double x, double **interpolationTable, int intepolationTblRows);
NUMERICLIBRARY_API enum DifferenceTableDirection { Forward, Backward, Center, Evaluate };
NUMERICLIBRARY_API DifferenceTableDirection ComputeDifferenceTable(double **interpolationTable, double **differenceTable, int tableRows, int tableCols, double x,
	DifferenceTableDirection direction = DifferenceTableDirection::Evaluate);
NUMERICLIBRARY_API double NewtonPolynomialInterpolation(double x, double **interpolationTable, int intepolationTblRows,
	DifferenceTableDirection direction = DifferenceTableDirection::Evaluate);
NUMERICLIBRARY_API double LagrangePolynomialInterpolation(double x, double **interpolationTable, int intepolationTblRows);

NUMERICLIBRARY_API double partialDerivative(double(*f)(double, double), double *varValues, int *changeIndex,
	int sizeOfchangeIndex, int orderOfDerivative);

NUMERICLIBRARY_API double EulersMethod(double x0, double y0, double x, double(*f)(double, double));
NUMERICLIBRARY_API double EulersMethodModified(double x0, double y0, double x, double(*f)(double, double)); NUMERICLIBRARY_API double EulersMethodModified(double x0, double y0, double x, double(*f)(double, double));
NUMERICLIBRARY_API double RungeKuttaFourthOrder(double x0, double y0, double x, double(*f)(double, double));
NUMERICLIBRARY_API void RungeKuttaForTwoDifferentialEquationsFourthOrder(double *x, double *y, double *z, int sizeOfVariableArray, double(*f[])(double, double, double), double h = .1);
NUMERICLIBRARY_API void AdamsBashforthMoultonFourthOrderForTwoDifferentialEquations(double x0, double y0, double z0, double x, double **xArr, double **yArr, double **zArr, int **sizeOfVariableArray, double(*f[])(double, double, double), double h);