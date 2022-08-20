

#include "physim.h"


using namespace physim; 


int main() {

    physics::position::coordinate coord(3.); 
    coord.print();
    std::cout << std::endl; 
    // physics::position::position pos(3., 4., 13); 
    // pos.print(); 

    // pos.convert_to(physics::units::defined::cm); 
    // pos.print(); 
    

    return 0;   
}