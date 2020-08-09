1) I should have used error codes as return values if I wanted to keep it to pure C.
But I used exceptions for error handling and that is a C++ concept.
But now i have manage the resource allocation and deallocation carefully as exceptions cause
the control to go out of scope without finishing the work being done.
http://stackoverflow.com/questions/385975/error-handling-in-c-code
https://msdn.microsoft.com/en-us/library/hh279678.aspx
2)I have also used C++ stl library vector in derivativeDividedByFactorialUsingSyntheticSubstitution.
I could have kept it to pure C by passing the array and its size but did not.
