

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.h"


using namespace physim; 
using namespace physics;
using namespace measurements;


int main() {

    measurement meas1(1, units::m); 
    measurement meas2(1, units::m); 
    measurement meas3(1, units::m); 
    measurement meas4(meas3 * 3); 
    meas4.print(); 
    std::cout << "\n"; 

    std::vector<measurement> pos1({meas1, meas2, meas4}); 
    for (auto i : pos1) i.print(); 
    
    std::vector<measurement> pos2(pos1 * 3); 
    for (auto i : pos2) i.print(); 

    // position p1(meas1, meas2, meas4);
    // p1.print_position();

    // position p2 = p1 * 3; 
    // p2.print_position(); 

    // assert(math::tools::are_close(p1.magnitude().value(), p1.distance(position(0, 0, 0)).value()));
    // p1.distance(position(0, 0, 0)).print();

    // assert(math::tools::are_close(p1.phi().value(), math::constants::pi / 4)); 
    // p1.phi().print();



    return 0; 
}