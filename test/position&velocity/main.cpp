

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the physics::position class
// last updated:    26/08/2022


#include "physim.hpp"

using namespace physics;
using namespace measurements;
using namespace units::SI;


int main() {

    position p1(measurement(2 * m), measurement(3 * m), measurement(4, m));
    p1.print_position();
    std::cout << "\n"; 

    position p2 = p1 * 3.; 
    p2.print_position(); 
    std::cout << "\n"; 

    position p3 = 3 * p1; 
    p3.print_position(); 
    std::cout << "\n"; 

    position p4 = p2 + p3; 
    p4.print_position(); 
    std::cout << "\n"; 

    position p5 = p4 / 3; 
    p5.print_position(); 
    std::cout << "\n"; 

    velocity v1(measurement(2 * mps), measurement(3 * mps), measurement(4, mps));
    v1.print_velocity();
    std::cout << "\n"; 

    velocity v2 = v1 * 3.; 
    v2.print_velocity(); 
    std::cout << "\n"; 

    velocity v3 = 3 * v1; 
    v3.print_velocity(); 
    std::cout << "\n"; 

    velocity v4 = v2 + v3; 
    v4.print_velocity(); 
    std::cout << "\n"; 

    velocity v5 = v4 / 3; 
    v5.print_velocity(); 
    std::cout << "\n"; 

    measurement meas(1, s); 
    meas.print_measurement();
    std::cout << "\n"; 

    position p6 = v5.get_velocity() * meas; 
    p6.print_position(); 
    std::cout << "\n"; 

    velocity v6 = p5.get_position() / meas;
    v6.print_velocity();
    std::cout << "\n"; 

    position p7 = p6 + v6.get_velocity() * meas;
    p7.print_position();
    std::cout << "\n"; 

    velocity v7 = v6 + p6.get_position() / meas;
    v7.print_velocity();
    std::cout << "\n"; 
    
    assert(math::tools::are_close(p1.magnitude().value(), p1.distance(position(0, 0, 0)).value()));
    p1.distance(position(0., 0., 0.)).print_measurement();
    std::cout << "\n"; 

    return 0; 
}