

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.hpp"


using namespace physics;
using namespace objects;


int main() {

    mass m1(30, units::kg); 
    m1.print_mass();  

    mass m2(30, units::kg);
    m2.print_mass();

    mass m3 = m1 + m2;
    m3.print_mass();

    mass m4 = m1 - m2;
    m4.print_mass();

    mass m5 = m3 / 3;
    m5.print_mass();

    mass m6 = m3 * 2; 
    m6.print_mass();

    charge c1(30, units::SI_derived::C); 
    c1.print_charge();  
    c1.print_position();
    c1.print_velocity();

    return 0; 
}