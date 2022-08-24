

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.h"


using namespace physim; 
using namespace physics;
using namespace objects;


int main() {

    mass m1(30, units::kg); 
    m1.print(); 

    mass m2(33.5); 
    m2.print(); 

    charge c1(30, units::SI_derived::C); 
    c1.print(); 

    charge c2(33.5); 
    c2.print(); 

    
    return 0; 
}