
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Test for physim/c++/include/math
// last updated:    03/07/2022


#include "../../include/math/integral.h"


int main() {


    Integral integral; 

    Cubic cub(3., 5., 0., 2.); 
    Sine sin(1., 6.); 

    Functor f('*', cub, sin); 

    Cosine cos(3., 0.5); 
    Logarithm log(M_E);

    Functor g('c', cos, log); 

    Functor k('/', f, g); 


    std::cout << "\n\nMidpoint" << std::endl; 
    integral.midpoint_fixed(3., 7., k, 1.e-10); 
    integral.print_integral(1.e-10);
    integral.print_error(); 
    
    std::cout << "\n\nTrapexoid" << std::endl; 
    integral.trapexoid_fixed(3., 7., k, 1.e-10); 
    integral.print_integral(1.e-10);
    integral.print_error(); 

    std::cout << "\n\nSimpson" << std::endl; 
    integral.simpson_fixed(3., 7., k, 1.e-10); 
    integral.print_integral(1.e-10);
    integral.print_error(); 
    
    // std::cout << "\n\nMean" << std::endl; 
    // integral.mean_fixed(3., 7., h, 1.e-7); 
    // integral.print_integral(1.e-7);
    // integral.print_error(); 
    
    
    return 0; 

}