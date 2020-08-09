#include "stdafx.h"
#include "NumericLibrary.h"

// This is the constructor of a class that has been exported.
// see NumericLibrary.h for the class definition
NotComputableException::NotComputableException(const char * message) : std::runtime_error(message)
{
}
