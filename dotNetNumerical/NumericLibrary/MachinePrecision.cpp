#include "stdafx.h"
#include "NumericLibrary.h"

NUMERICLIBRARY_API void GetMachinePrecision() {
	printf("%30s: %g\n", "FLT_EPSILON", FLT_EPSILON);
	printf("%30s: %g\n", "FLT_MIN", FLT_MIN);
	printf("%30s: %g\n", "nextafterf(0.0, 1.0)", nextafterf(0.0, 1.0));
	printf("%30s: %g\n", "nextafterf(1.0, 2.0)-1", (nextafterf(1.0, 2.0) - 1.0f));
	puts("");
	printf("%30s: %g\n", "DBL_EPSILON", DBL_EPSILON);
	printf("%30s: %g\n", "DBL_MIN", DBL_MIN);
	printf("%30s: %g\n", "nextafter(0.0, 1.0)", nextafter(0.0, 1.0));
	printf("%30s: %g\n", "nextafter(1.0, 2.0)-1", (nextafter(1.0, 2.0) - 1.0));
	puts("");
	printf("%30s: %Lg\n", "LDBL_EPSILON", LDBL_EPSILON);
	printf("%30s: %Lg\n", "LDBL_MIN", LDBL_MIN);
	printf("%30s: %Lg\n", "nextafterl(0.0, 1.0)", nextafterl(0.0, 1.0));
	printf("%30s: %Lg\n", "nextafterl(1.0, 2.0)-1", (nextafterl(1.0, 2.0) - 1.0));
}