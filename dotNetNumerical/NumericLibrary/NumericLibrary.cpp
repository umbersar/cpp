// NumericLibrary.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "NumericLibrary.h"


// This is an example of an exported variable
NUMERICLIBRARY_API int nNumericLibrary=0;

// This is an example of an exported function.
NUMERICLIBRARY_API int fnNumericLibrary(void)
{
    return 42;
}

// This is the constructor of a class that has been exported.
// see NumericLibrary.h for the class definition
CNumericLibrary::CNumericLibrary()
{
    return;
}
