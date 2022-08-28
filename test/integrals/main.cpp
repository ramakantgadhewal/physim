

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.hpp"
#include "gplot++.h"


using namespace physim; 


int main() {

    math::tools::integral integral; 
    math::functions::sine sin;

    integral.midpoint(0., math::constants::pi, sin, 10);
    assert(utilities::tools::are_close(integral.value(), 2.0082484079079745, 1.e-16));
            
    integral.midpoint(0., math::constants::pi, sin, 100);
    assert(utilities::tools::are_close(integral.value(), 2.000082249070986, 1.e-15));
            
    integral.midpoint(math::constants::pi, 0., sin, 10);
    assert(utilities::tools::are_close(integral.value(), -2.0082484079079745, 1.e-16));
    
    integral.midpoint(0., 1., sin, 10);
    assert(utilities::tools::are_close(integral.value(), 0.45988929071851814, 1.e-17));
    
    integral.midpoint(1., 2., sin, 30);
    assert(utilities::tools::are_close(integral.value(), 0.9564934239032155, 1.e-15));
     
    std::cout << "midpoint test passed successfully \n"; 


    integral.simpson(0., math::constants::pi, sin, 10); 
    assert(utilities::tools::are_close(integral.value(), 2.0001095173150043, 1.e-16));  

    integral.simpson(0., math::constants::pi, sin, 100); 
    assert(utilities::tools::are_close(integral.value(), 2.000000010824504, 1.e-15));  

    integral.simpson(0., 1., sin, 10); 
    assert(utilities::tools::are_close(integral.value(), 0.45969794982382056, 1.e-17));  

    integral.simpson(1., 2., sin, 30); 
    assert(utilities::tools::are_close(integral.value(), 0.9564491489761575, 1.e-16));  
    
    std::cout << "simpson test passed successfully \n"; 


    integral.trapexoid(0., math::constants::pi, sin, 10); 
    assert(utilities::tools::are_close(integral.value(), 1.9835235375094546, 1.e-16));  

    integral.trapexoid(0., math::constants::pi, sin, 100); 
    assert(utilities::tools::are_close(integral.value(), 1.9998355038874436, 1.e-16));  

    integral.trapexoid(0., 1., sin, 10); 
    assert(utilities::tools::are_close(integral.value(), 0.45931454885797635, 1.e-17));  

    integral.trapexoid(1., 2., sin, 30); 
    assert(utilities::tools::are_close(integral.value(), 0.956360580669458, 1.e-15));  

    std::cout << "trapexoid test passed successfully \n"; 



    return 0;
}