

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.hpp"


using namespace physim; 
using namespace physics;
using namespace measurements; 

int main() {


    measurement meas1(34, units::SI_derived::Hz); 
    meas1.print(); 

    return 0; 
}