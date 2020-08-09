#include "stdafx.h"
#include "NumericLibrary.h"

//PDE are divided into elliptic, parabolic and hperbolic and can be first order or higher order.

//For PDE, the example show below is a bounday value problem(BVP) for elliptic PDE. Just like ODE BVP, PDE BVP can also makes use of finite 
//difference approximationss for partial derivatives and then computes the solution at mesh points. And then we just improve the
//solution by using a iterative method. 
//But in ODE, after using finite difference approximations to the derivatives occuring in the eq, we solve the resulting system of linear
//equations using any standard procedure (like LUDecompositionLinearEquationsSolver coded in MatrixLinearAlgebra.cpp)
//
//what about initial value problems for PDE

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

//Lets start with a example. We have a laplace equation(elliptic PDE):Uₓₓ + Uᵧᵧ = 0.  Use the finite difference approximations as shown in
//PartialDerivatives.cpp: 
//(fᵢ₊₁,j - 2*fᵢ,j + fᵢ₋₁,j) / 2*h +  (fᵢ,j₊₁ - 2*fᵢ,j + fᵢ,j₋₁) / 2*k = 0. Assume h=k
//(fᵢ₊₁,j - 2*fᵢ,j + fᵢ₋₁,j) / 2*h +  (fᵢ,j₊₁ - 2*fᵢ,j + fᵢ,j₋₁) / 2*h = 0
//((fᵢ₊₁,j - 2*fᵢ,j + fᵢ₋₁,j) +  (fᵢ,j₊₁ - 2*fᵢ,j + fᵢ,j₋₁)) / 2*h = 0
//((fᵢ₊₁,j - 2*fᵢ,j + fᵢ₋₁,j) +  (fᵢ,j₊₁ - 2*fᵢ,j + fᵢ,j₋₁)) =0
//fᵢ₊₁,j + fᵢ₋₁,j + fᵢ,j₊₁ + fᵢ,j₋₁ - 4*fᵢ,j = 0
//fᵢ,j = (fᵢ₊₁,j + fᵢ₋₁,j + fᵢ,j₊₁ + fᵢ,j₋₁ )*(1/4). This is called the standard 5 point formula
//Instead of the standard 5 point formula, we can also use diagonal 5 point formula:
//fᵢ,j = (fᵢ₋₁,j₋₁ + fᵢ₊₁,j₋₁ + fᵢ₊₁,j₊₁ + fᵢ₋₁,j₊₁ )*(1/4).

//now we wish to solve Uₓₓ + Uᵧᵧ = 0 in a bounded region R with boundary C and the value of f is specified everywhere on C. we assume
//the h=k.                    
//Let c1=0, c2=500, c3=1000, c4=500, c5=0, c6=1000 and so on as shown below

//                       500       1000      500
//             c13       c12       c11       c10          
//       0      ●---------●---------●--------●---------●c9 0              
//              |         |         |         |         |       
//              |         |         |         |         |       
//              |         |f₇       |f₈       |f₉       |       
//    1000  c14 ●---------●---------●---------●---------●c8 1000              
//              |         |         |         |         |       
//              |         |         |         |         |       
//              |         |f₄       |f₅       |f₆       |       
//    2000  c15 ●---------●---------●---------●---------●c7 2000              
//              |         |         |         |         |       
//              |         |         |         |         |       
//              |         |f₁       |f₂       |f₃       |       
//    1000  c16 ●---------●---------●---------●--------●c6 1000              
//              |         |         |         |         |       
//              |         |         |         |         |       
//              |         |         |         |         |       
//              ●---------●---------●---------●--------●              
//              c1        c2        c3        c4        c5 
//              0        500       1000      500        0
//          
//          
//Now since the values we have are the values on the boundary of the mesh, so we have to make use of the boundary points first to 
//calculate the unknowns. Now we do not have to restrict the application of the application of the 5 point formulas to immediate 
//neighbours. So, f₅ is:
//f₅=(2000+2000+1000+1000)*1/4     =1500   standard 5 point
//f₇=(0+1500+1000+2000)*1/4        =1125   diagonal 5 point. we use the value of f₅ computed above
//f₉=(0+1500+1000+2000)*1/4        =1125   diagonal 5 point. we use the value of f₅ computed above
//f₈=(1000+1500+1125+1125)*1/4     =1188    
//f₁=(0+1500+1000+2000)*1/4        =1125   diagonal 5 point. we use the value of f₅ computed above
//f₄=(1125+2000+1125+1500)*1/4     =1438   standard 5 point  
//f₃=(0+2000+1500+1000)*1/4        =1125   diagonal 5 point. we use the value of f₅ computed above
//f₆=(1125+2000+1500+1125)*1/4     =1438   standard 5 point  
//f₂=(1125+1125+1500+1000)*1/4     =1188   standard 5 point   
//
//so now we have all the values at mesh points.note that due to symmetry in the problem f₄=f₆, f₁=f₃, f₇=f₉, f₇=f₁, f₈=f₂ and f₉=f₃.
//So we only need to evaluate f₇, f₈, f₄ and f₅.
//To recapture, what we did was that we replaced the partial derivatives by finite difference approximations and then obtainied 
//the function values at mesh points. Once we have the values at the mesh points, we improve their accuracy by iterative methods. 
//iterative methods are basically the standard 5 point formula. Gauss sidel iterative formulas are(they use the latest iterative 
//value if available):
//fⁱ⁺¹₇=(1000+500+fⁱ₈+fⁱ₄)*1/4          standard 5 point.
//fⁱ⁺¹₈=(1000+fⁱ⁺¹₇+fⁱ⁺¹₇+fⁱ₅)*1/4      standard 5 point.
//fⁱ⁺¹₄=(2000+fⁱ⁺¹₇+fⁱ⁺¹₇+fⁱ₅)*1/4      standard 5 point.
//fⁱ⁺¹₅=(fⁱ⁺¹₈+fⁱ⁺¹₈+fⁱ⁺¹₄+fⁱ⁺¹₄)*1/4   standard 5 point.

//First iteration:i=0
//f¹₇=(1000+500+1188+1438)*1/4		=1032       standard 5 point.
//f¹₈=(1000+1032+1032+1500)*1/4     =1141       standard 5 point.
//f¹₄=(2000+1032+1032+1500)*1/4     =1391       standard 5 point.
//f¹₅=(1141+1141+1391+1391)*1/4     =1266       standard 5 point.
//
//