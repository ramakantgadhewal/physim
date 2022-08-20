

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.h"
#include "gplot++.h"


using namespace physim; 
using namespace physics;


int main() {

    measurements::measurement meas1(33, units::SI::m); 
    measurements::measurement meas2(3, units::SI::s); 
    measurements::measurement meas3 = meas1 / meas2; 

    meas1.print(); 
    meas2.print(); 
    meas3.print(); 


    return 0; 
}