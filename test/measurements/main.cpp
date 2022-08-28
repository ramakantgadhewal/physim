

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.hpp"


using namespace physim; 
using namespace physics;
using namespace measurements;


int main() {

    std::cout << "\nmeasurements\n";
    measurement meas1(336, units::m); 
    measurement meas2(3, units::SI_derived::mps); 
    measurement meas3 = meas1 / meas2; 

    meas1.print(); 
    meas2.print(); 
    meas3.print(); 

    std::cout << "\nfixed_measurements\n";
    fixed_measurement fmeas1(24, units::km); 
    fixed_measurement fmeas2(340, units::m); 
    fixed_measurement fmeas3 = fmeas1 + fmeas2; 
    fixed_measurement fmeas4 = fmeas2 + fmeas1; 

    fmeas1.print(); 
    fmeas2.print(); 
    fmeas3.print(); 
    fmeas4.print(); 

    std::cout << "\nuncertain_measurements\n";
    uncertain_measurement umeas1(15.6, 0.3, units::m.pow(2));
    uncertain_measurement umeas2(23.7, 0.6, units::m);
    uncertain_measurement umeas3 = umeas1 * umeas2;
    
    umeas1.print();
    umeas2.print();
    umeas3.print();
    std::cout << "\n";

    return 0; 
}