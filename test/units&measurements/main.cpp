

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.hpp"


using namespace physics;
using namespace measurements;
using namespace units; 


int main() {

    std::cout << "\nmeasurements\n";
    measurement meas1(336 * m); 
    measurement meas2(3, SI_derived::mps); 
    measurement meas3 = meas1 / meas2; 

    meas1.print(); 
    meas2.print(); 
    meas3.print(); 

    std::cout << "\nuncertain_measurements\n";
    uncertain_measurement umeas1(15.6, 0.3, m.pow(2));
    uncertain_measurement umeas2(23.7, 0.6, m);
    uncertain_measurement umeas3 = umeas1 * umeas2;
    
    umeas1.print();
    umeas2.print();
    umeas3.print();
    std::cout << "\n";

    return 0; 
}