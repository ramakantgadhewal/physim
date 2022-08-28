
#include "../../physim.hpp"


int main() {

    utilities::timer timer; 
    timer.start(); 
    std::vector<double> v = utilities::read_from_file<double>("data.dat"); 

    timer.start(); 
    assert(utilities::tools::are_close(math::descriptive_statistics::mean(v), 30.23231, 1.e-5)); 
    timer.pause(); 
    timer.print();

    timer.start();
    assert(utilities::tools::are_close(math::descriptive_statistics::variance(v), 282326.76577, 1.e-5));
    timer.pause(); 
    timer.print();

    timer.start();
    assert(utilities::tools::are_close(math::descriptive_statistics::sd(v), 531.34430, 1.e-5));
    timer.pause(); 
    timer.print();

    timer.start();
    assert(utilities::tools::are_close(math::descriptive_statistics::median(v), 12.74255, 1.e-5));
    timer.pause(); 
    timer.print();

    return 0; 
}