#include "stdafx.h"

//Comment the file to compile project

//1) try to always use explicit with constructors(with parameters) to prevent unexpected and unintended implicit type conversions.
class A {
public:
	explicit A(int x);
};

//2) copy constructor is used to initialize an object with a different object of the same type and copy assignment operator is used to copy the
//value from one object to another object of the same type.
class widget {
public:
	widget();//default constructor
	widget(const widget& rhs);//copy constructor
	widget& operator=(const widget& rhs);//copy assignment operator
};

widget w1; //invoke default constructor.. if we had written widget w1();, then it would have meant we are declaring a function w()
widget w2(w1);//invoke copy constructor
widget w3 = w1;//invoke copy constructor
//copy constructor is also called when an object is passed-by-value(instead of by reference) to a function.

w1 = w2;//invoke copy assignment operator

//3) underfined behaviour
int *p = 0;// null pointer
std::cout << *p;//underfined behaviour

//4) there is no equivalent of C# interface in c++

//5)prefer const, enum and inline to #define aka prefer compiler to pre processor. 
#define ASP_RATIO 53.1
const double Asp_ratio = 53.1;//replaced the macro with const
//remember in case  of pointers, not only the pointer be a constant but what it points to be as well
const char* const myName = "troll";//what is difference between this and const char const *myName;

char * const a; //*a is writable, but a is not; in other words, you can modify the value pointed to by a, but you cannot modify a itself.  a is a constant pointer to char.
const char * a; //a is writable, but *a is not; in other words, you can modify a (pointing it to a new location), but you cannot modify the value pointed to by a. same as below.
char const * a;//similar as above, a is a pointer to a const char.

//class level constants are better declared static as well
class widget {
private:
	static const int Num = 3;
};

//in addition to defining global constant values, define macros are also used to provide function like behaviour without the actual cost of function call. In those cases use inline funtions!

//6) use const whenever possible. see above for const pointers. iterators aree modelled on pointers. Declaring a iterator const is like declaring a pointer const, i.e., 
//   T * const ptr. So the iterator isn't allowed to point to something different but the thing it points to may be modified. if you want a iterator that points to something that can't be 
//   modified(const T * ptr), then use const_iterator
std::vector<int> vec;

const std::vector<int>::iterator it = vec.begin(); // like T * const
*it = 10; // ok
++it; //error. it is const

std::vector<int>::const_iterator cit = vec.begin(); // like const T * 
*cit = 10; // error
++cit; //ok

//use const whenever possible for local variables, parameters, function return value and for member functions itself. const member functions do not alter the state of the object. Thus
//const objects(as well as non-const objects I think) can call the const member functions. Member function differing in only const'ness can be overloaded.
class TextBlock {
private:
	std::string data;
public:
	const char& operator[](const std::string::size_type pos) const { return data[pos]; }//operator for const objects
	char& operator[](const std::string::size_type pos) { return data[pos]; }//operator for non-const objects
};

TextBlock t("some randon things");
std::cout << t[1];//uses non-const operator[] func

const TextBlock ct("some randon things");
std::cout << ct[1];//uses  const operator[] func

//mutable data members are those which can be modified in const member functions(so const member funcs are changing the state of the object)

//7) try to always initialize built in types either when defining or by reading value into them from a stream.
// and for custom types, make sure constructor initialize everything in the class. But do not confuse initialization with assignment. Assignment happens in the
// constructor body whereas initialization happens before it
class ex {
private:
	int i;
	std::vector<int> a;
	std::string j;
public:
	ex() :i(0), a(), j() {

	}
	ex(int q, const std::vector<int>& w, const std::string& e) :i(q), a(w), j(e) {

	}
};

//before the constructor body in entered, all the custom type data members have their constructors executed either manually as we are doing OR if we had not done so, their default constructors would be called
// automatically. For built-in type int, if we forgot to initialize it, there is no guarantee it would have been initialized automatically.
//base class constructor is called before derived class cnstr. the order in which the data members are initialized(either through manual initialization list or through automatic call to default constructors)
// is the the order in which the members are declared in the class. so to avoid confusion, always initialize members in the same order in which they are
//declared.

//8) compiler will generate a default constructor, copy constructor, assignment operator and destructor if we do not declare them but try to use them.
class def {

};
//will generate something like
class def {
	def() { .. }
	def(const def& rhs) { ... }
	~def() { ... }
	def& operator=(const def& rhs) { ... }
};
//the compiler generate copy cnstr and assignment operator to do a simple copy the data members of source object to target object. So if your class
// has pointers, the address is copied. so both the copied oject and source object are pointing to the same object through the pointer(shallow copy)
//where as what should have been done is that it should have cloned the object pointer to by pointer and use it in the copied object.

//also if you have only declared a cnstr with arguments, the compiler won't generate a default cnstr.

//9) For some classes, you do not want to define copy cnstr and assignment operator as you want the objects to be unique. In that case declare(but do not define them) private versions of the functions
// so that clients cant call them. In that case compilers also wont implement them as they are already declared in the class.
class unq {
private:
	unq(const unq&);
	unq& operator=(const unq&);
private:
	int *ptr;
};
unq u1;
unq u2(u1);//compile time or linker error depending on how you call it.
u1 = u2;//compile time or linker error depending on how you call it.

//10) if you are using virtual funcs, then after you are done with the virtual func calls, you have to delete the pointer to release
//the memory after you are done with. Now if you delete the ptr, it only delete the base class object. So you also have to make the destructor virtual so as to call the correct destructor:
class shape {
public:
	shape();
	virtual ~shape();
};
Shape *ptr;
ptr->draw();
delete ptr;//unless destructor is declared virual in base class, this call will not execute the destructor of the derived class.

//so if a class has virtual funcs, then declare the destructor as virtual

//also note that ptr can be set using 2 techniques:
ptr = new Rectangle();//create a new anonymous object Rectangle and set the address of it in ptr
//OR if the object already exists, then by just taking it's address
ptr = &rec;


//11) prevent exceptions from leaving destructors. So either swallow the exception in the destructor or abort the program in catch block . and if the client needs to be given
//the oppurtunity to be able to react to the exceptions, then the move the desctruction logic to a callable function(that client can call) and call that func from within
//destructor as a backup otion

//12) do not call virtual functions in cnstrs and dstrs as they would not be forwarding the call to the more derived classes

//13) have assignment operators return a reference to *this (because assignment yields the lhs).

//14) handle assignment to self in operator=. Techniques used include comparing addresses of source and target objects, careful statement ordering and copy and swap.
unq::unq& operator=(const unq& rhs) {
	if (this != &rhs) {
		delete ptr;
		ptr = rhs.ptr;
	}
	return *this;
}

unq::unq& operator=(const unq& rhs) {
	int *ptrBackup = ptr;
	ptr = rhs.ptr;
	delete ptrBackup;

	return *this;
}


//15) copying funcs(copy cnstr and assignment op.) should copy all the data members(including those of base classes in case of inheritance hierarchy)

//15.1) do not delete a pointer more than once as it is undefined behaviour.

//16) use objects to manage resources. Mean client code can't be trusted with deletion of dynimcally allocated memory. wrap the ptr in a class so that the destructor
// deletes the ptr. U can use the handle class approach as shown in Accelerated C++ book. Some newer approaches are to use auto_ptr(smart pointer). using auto_ptr, we wrap the 
// raw pointer in a smart pointer class whose destructor calls the delete on raw pointer. So as soon as smart ptr goes out of scope, the memory is freed.

//client managing the deletion
void f() {
	shape *ptr = CreateShapeObjects();
	...
		delete ptr;
}

//auto_ptr smart pointer managing the deletion. auto_ptr is deprecated. use unique_ptr instead.
void f() {
	std::auto_ptr<shape> ptr(CreateShapeObjects());
	...
}

//there is one thing to remember about auto_ptr. It is that if we copy construct or assign a auto_ptr to another, then the source auto_ptr is set to null(null means 0..note that it is only set to null, not deleted)
//In the new standard auto_ptr is replaced by unique_ptr. Also since

//shared_ptr smart pointer managing the deletion. But with shared_ptr we will get the expected copy and assignment behaviour. It is a reference couting
//shared pointer. So it will keep track of how many references to the underlying resource exists and automaticall deletes when nobody is pointing to it
//so it is almost like garbage collection but with a glitch that it can break circular references and that can cause memory leak.
//use std::make_shared or std::make_unique whereever possible instead of new'ing an object(you can use them when using custom deleters
//refer item 16 down below)
void f() {
	std::shared_ptr<shape> ptr(CreateShapeObjects());
	...
		std::shared_ptr<shape> ptr2(ptr);
	ptr = ptr2;
	...
}

//keep in mind that smart pointers are just like the handle classes used in Accelerated C++ book. So smart pointers also overload the -> and * operators.
void f() {
	std::shared_ptr<shape> ptr(CreateShapeObjects());

	ptr->Draw();
	(*ptr).Draw();
	...
}

//17) Keep the same form when using new and delete. If allocating an array with new, do not forget to use delete[].
std::string *s = new std::string[10];
delete[]s;

//remember that new allocates memory and then call constructor(s) of object(s) for which memory was allocated. So it initializes memory as well something 
//that malloc does not do. Now deletes calls destructors for object(s) and then makes the memory avaiable to the system for resuse.
//funny pattern: do your deletes in destructors -- delete calls destructors AND do you new in constructors -- new calls constructors 

//18)make interface hard to use incorrectly:
class Date {
public:
	Date(int day, int month, int year);
};

Date dt(12, 27, 1999);

struct Day
{
	Day(int d) :val(d) {}
	int val;
};

struct Month
{
	Month(int m) :val(m) {}
	int val;
};

struct Year
{
	Year(int m) :val(m) {}
	int val;
};
class Date {
public:
	Date(const Day& d, const Month& month, const Year& year);
};

Date dt(Day(1), Month(1), Year(1999));

//now we can also restrict the values that user provides. one way could be to use enum or in the setter functions we can validate the values passed

//19)prefer pass-by-reference-to-const TO pass-by-value for non native types(class, structs) but for native types (int, double, char), iterators and function objects prefer to pass-by-value. 
//Also remember that references are actually internally implemented as pointers. So the virtual functions polymorphic behaviour can be exploited both by pointers 
//as well references.

//pass-by-reference-to-const for class objects prevents the cost of constructing a copy of the object(through copy constructors) and they save a bit of memory as well. 
//But it is more to do with construction cost

void ex(int q, const std::vector<int>& w, const std::string& e) {

}

//20) don't try to return a reference when you must return an object. That is because the local objects would be destroyed once the scope of function finishes
// you return a reference in the assignment operator because the object being returned is 'this' which exists outside the scope of the assignment function.

//local variables are created on stack. Heap based objects are created through use of 'new'

//21) in C# u declare utlitiy funcs as static members of c utlity class. We can do that in c++ as well but it is more natural for them to be just in the same namespace(but 
//maybe declare in separate header file), not in any utility class

//22) minimize casting: 
//c style casts:
(T)expression;//cast expression to be of type T

//new style casts. This needs review as I do not think it is correct:
const_cast<T>(expression);//cast away const'ness of objects
dynamic_cast<T>(expression);//primarily used to perform safe downcasting 
reinterpret_cast<T>(expression);//for reinterpreting an object as something else
static_cast<T>(expression);//from non-const to const, int to double, void* to typed pointers...so in a way it is doing sort of 'upcasting'

int a = 7;
double* p1 = (double*)&a;			// ok (but a is not a double)
double* p2 = static_cast<double*>(&a);	// error
double* p2 = reinterpret_cast<double*>(&a);	// ok: I really mean it

const int c = 7;
int* q1 = &c;			// error
int* q2 = (int*)&c;		// ok (but *q2=2; is still invalid code and may fail)
int* q3 = static_cast<int*>(&c);	// error: static_cast doesn't cast away const
int* q4 = const_cast<int*>(&c);	// I really mean it

int x, y;
double d = static_cast<double>(x) / y;//c++ new style cast
double d1 = (double)x / y;//c style cast

//23) defining a function inside the class definition is a request to inline the function whereever it is called. You can define the function outside
//and append 'inline' to the function definition. 
//inline functions and templates are generally defined in header files.

//24) public inheritance means "is-a" relationship. The other relationships that can exist between classes are "has-a" and "is-implemented-in-terms-of"

//25) avoid hiding inherited names. one example, do nt redefine non-virtual funs in derived classes as they hide the base funcs.
// if you have to, then the func has to be virtual(not a non virtual).

//26) Do not redefine a virtual functions's inherited default parameter value (we are only talking about virtual func redefinition
// as we should not be re-defining non-virtual funcs..see above). default params are statically bound(bound at compile time, not run time). So the base class default params would be used
//even if the derived class virtual func is being called.

//27)public iheritance means 'is-a'. Composition means 'has-a' or 'is-implemented-in-terms-of' 
//composition is the realtion which arises when one type contains objects of another type. So when a class Person defines a sting member name, it is composition. Person 'has-a' name.

//as a guide, public inheritance or 'is-a' suggests that all the properties of base class are valid for derived class as well
//so you should be able to pass in derived object to base class variables to get 'is-a' relation to succeed.
//think about it. If you want to construct class Set using List, then instead of public inheritance('is-a'), we should use composition and enacpsulate a
//List object inside of Set because not all properties of List stand truu for Set. example: items can be duplicated in List but not in Set.

//also square should be composed of rectangle, instead of being derived from rectangle. 

//a good test is that all the data members of base should be valid data members of derived. And all the member funcs of base should be valid for derived as well
//now the reason rectangle cant be base class of shape is that it has a func 'makebigger' which takes 2 params length and breadth.
//now we can try to make this func virtual as well, but still square version of func will still require 2 params which does not make sense.
//now we could add some validation in square's makebigger to impose that length should be equal to breadth and that could work. But then you should remember to pass
//references or pointers around for square objects. If you set square object to a rectangle variable(local or a function param) which is not a reference or pointer,
//then it might call the func makebigger and that will call rectangle's version to get around the restrictions set on square.

//28) private inheritance also means 'is-implemented-in-terms-of' just like composition. also a privately inherited class object can't be
//set to a base class variable. 
//Public and Protected members inherited from base class become private members in the derived class. So derived class objects
//can't even call the public members from base class(can only call from inside the derived class).

//now if both composition and private inheritance mean 'is-implemented-in-terms-of' , how do we decide to use which?
//use composition whenever you can and private inheritance whenever you must.

//29) there is another type of inheritance, virtual. WTF!

//30) try to avoid multiple inheritance(MI). If not possible, then implement it as C# does. Use interfaces for MI. Now C++ does not have interfaces.
//so use abstract classes(But should it have data members? Are data members allowed in C# interfaces).

//31)for templates, template<class T> is interchangeable with template<typename T>. when we have to identify nested dependent type name(dependent on T)
//(or maybe even refer to a static var inside T?), we have to use 'typename' prefix. But 'typename' cant be used inside base lists or member initialization lists

//32)template specialization: when the template class or template function needs to behave differently for a particular type.

//33) template classes are instansiated during compile time. Factor template parameter independent code out of templates to avoid code bloat.

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =

//Effective Modern C++

//1) parameters are lvalues bu the arguments with which they are initialized may be rvalues or lvalues.
//another way to remember is that lvalue is a value which has a name(variable name) but rvalue is not named

void someFunc(widget w);
widget wid;
someFunc(wid);//in this call to func w is a copy of wid that is created via copy construction
someFunc(std::move(wid));//in this call to func w is a copy of wid that is created via move construction

//tatt:old c++ style references are actually lvalue references(formed by placing an & after some type). new rvalue reference is formed by 
//placing an && after some type
A a;
A& a_ref1 = a;  // an lvalue reference
A&& a_ref2 = a;  // an rvalue reference

//An rvalue reference behaves just like an lvalue reference except that it can bind to a temporary(an rvalue), whereas you can not bind a(non const) lvalue reference to an rvalue.
A&  a_ref3 = A();  // Error!
A&& a_ref4 = A();  // Ok

//in some contexts && means rvalue references and in some contexts it means universal references.

//2) function objects created through lambda expressions are known as closures.

//3) universal references:
void someFunc(widget&& w); //w is a universal reference

//4) decltype. given a name or an expression, decltype tells you it's type
const int i = 0;
decltype(i);// const int

//5) of c++'s three ways to do initialization(=,{},()), only the braces work in all situations. So they are called 
//uniform initialization
widget cp;
widget w1(cp);
widget w2 = cp;
widget w3{ cp };

widget w4();//this is taken to be a function declaration. so wrong
widget w5;//default constructor
widget w6{};//default constructor

std::vector<int> sas1(2, 1); //it creates a 2-element vector with initial values 1 for both
std::vector<int> sas2{ 2,1 };//it creates a 2-element vector with initial values 2 and 1

//there are some caveats with {}. So I would keep it simple. For initializing simple types(int,double), I would use = for initialization
//for custom types, I would use () for initilization unless I have to initialize a vector where I would use {}.

//6) prefer nullptr to 0 and NULL

//7)prefer alias over typedef for defining synonyms
typedef std::unique_ptr<std::map<std::string, std::string>> shortName;//typedef to reduce the amount of typing
using shortName = std::unique_ptr<std::map<std::string, std::string>>;//same thing using alias

typedef void(*fp)(int, const std::string&);//typedef for a function pointer takes int and string and returns void
using fp = void(*)(int, const std::string&);//same using alias declaration

//7)use scoped enums. enums leak. Means although names declared inside enum clurly braces seem to be scoped, they are not
enum Color { White, Black, Red };//White is available outside braces as well. So can't declare variable white!
bool White = false;//won't compile as white already declared in this scope

enum class Color { White, Black, Red };//scoped enum...now u can declare White variable
//also with scoped enums, you won't be able to use enum values interchangeably with int values(for example in comparisons like White>3)
//you would have to use casting to do same. if(static_cast<double>(Color::White)>3)...

//8) Instead of making some functions unavailable by making making them private(with just a declaration but no definition),
//use "= delete" for same functionality

//old way of making the assignment operator and copy constructor unavaiable to prevent copying and assignment.
//We had to declare them private. If we had not declared them, they would have been added automatically 
//by the compiler
class widgetOld {
public:
	widget();//default constructor
private:
	widget(const widget& rhs);//make copy constructor unavailable
	widget& operator=(const widget& rhs);//make copy assignment operator unavailable
};

//modern way of preventing copying and assignment

class widgetNew {
public:
	widget();//default constructor
	widget(const widget& rhs) = delete;//make copy constructor unavailable
	widget& operator=(const widget& rhs) = delete;//make copy assignment operator unavailable
};

//9) delare overriding functions(exploiting virtual behaviour) with override to prevent overloading by mistake
class base {
	virtual void mf();
	virtual void mf1() const;
	virtual void mf2(int i);
};

class derived :public base {
	virtual void mf();//this is overriding correctly so you can exploit he virtual behaviour correctly
	virtual void mf1();//this is not overriding but overloading. so you cannot exploit he virtual behaviour correctly
	virtual void mf2(unsigned int i);//this is not overriding but overloading. so you cannot exploit he virtual behaviour correctly
};

//the above code will compile but had my intention was to override and thus it is a bug which can be prevented by 
//declaring the funcs as 'override' 
class derivedImproved :public base {
	virtual void mf() override;//this is overriding correctly so you can exploit he virtual behaviour correctly
	virtual void mf1()override;//this line won't compile now
	virtual void mf2(unsigned int i)override;//this line won't compile now
};

//10) final applied to virtual funcs prevents them from being overridden in derived classes. And if final is applied
//to a class, the class can't be used as a base class

//11) prefer const_iterators ot iterators whenever you do not need to modify what iterator points to

//12) declare functions noexcept if they won't emit exceptions

//13) mutable types can be used in const functions???

//14) std::atomic and std::mutex for concurrent scenarios.

//15) point number 8 from old c++(compiler will generate a default constructor, copy constructor, assignment operator and destructor if we do not declare them but try to use them.)
//C++ 11 adds two more, move constructor and move assignment operator. Again, they would only be generated only if we do not declare them but try to use them.
class def {

};
//will generate something like
class def {
	def() { .. }
	def(const def& rhs) { ... }
	~def() { ... }
	def& operator=(const def& rhs) { ... }

	def(def&& rhs) { ... };//move constructor
	def& operator=(def&& rhs) { ... };//move assignment operator
};

//16) both shared_ptr and unique_ptr support custom deleters...So for example, you can use custom deleters for logging
//whenever an delete is called(internally both these smart pointers do call delete for resource deallocation)

//17)Lambda expressions is just an expression
std::find_if(sas1.begin(), sas1.end(), [](int val) {return 0 < val && val < 10; });

//closure is a object created by lambda. Closure holds copies of or references to captured data. the 3rd argument above is an closure
int x;
auto c1 = [x](int y) {return x*y > 55; };
auto c2 = c1;//copy of c1

//a by reference capture causes a closure to contain a reference to a local variable or to a parameter that's available in
//the scope where the lambda is defined. if the lifetime of a closure created from that lambda exceeds the lifetime of the 
//local variable or parameter, the reference in the closure will dangle

auto c1 = [&x](int y) {return x*y > 55; };

//17) prefer task based programming to thread based
int doAsynWorkFunc();

std::thread t(doAsynWorkFunc);//thread based async call...it offers noway to get return values from asynchronously run funcs
//and u have manually mange thread oversusbsription, load balancing

auto fut = std::async(doAsynWorkFunc);//task based async call