1) I should have used error codes as return values if I wanted to keep it to pure C.
But I used exceptions for error handling and that is a C++ concept.
But now i have manage the resource allocation and deallocation carefully as exceptions cause
the control to go out of scope without finishing the work being done.
http://stackoverflow.com/questions/385975/error-handling-in-c-code
https://msdn.microsoft.com/en-us/library/hh279678.aspx
2)I have also used C++ stl library vector in derivativeDividedByFactorialUsingSyntheticSubstitution.
I could have kept it to pure C by passing the array and its size but did not.
3) I have used function overloading which is not part of C standard.



Example codes: 
Python http://www.cs.gordon.edu/courses/mat342/python.html
Matlab https://www.youtube.com/channel/UCq41scHkyvxM5Ku6hJAkHGg Montana state univ
Matlab https://math.berkeley.edu/~mgu/MA128A2008F/
c      http://sepwww.stanford.edu/sep/sergey/128A/


TODO: Monte carlo methods and fourier transform
