

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    20/08/2022


#include "physim.h"


using namespace physim; 
using namespace physics;
using namespace measurements;
using namespace tools; 


int main() {

    measurement meas1(1, units::m); 
    measurement meas2(1, units::m); 
    measurement meas3(1, units::m); 

    position pos(meas1, meas2, meas3);
    pos.x().print(); 
    pos.print();

    assert(math::tools::are_close(pos.magnitude().value(), pos.distance(position(0, 0, 0)).value()));
    pos.distance(position(0, 0, 0)).print();

    assert(math::tools::are_close(pos.phi().value(), math::constants::pi / 4)); 
    pos.phi().print();

    return 0; 
}