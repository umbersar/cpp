#include "stdafx.h"
#include "NumericLibrary.h"

// This is the constructor of a class that has been exported.
// see NotImplementedException.h for the class definition
NotImplementedException::NotImplementedException(const char * message = "Functionality not yet implemented!") : std::logic_error(message)
{
}
