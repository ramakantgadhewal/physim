
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Testing Mass(class) and its gravitational functions.
// last updated:    11/07/2022


#include "../../include/physics/tools/mass.h"
#include "../../include/physics/tools/coordinates.h"


int main() {

    Mass m(100, "k"); 
    m.print_mass();  
    
    return 0; 
}