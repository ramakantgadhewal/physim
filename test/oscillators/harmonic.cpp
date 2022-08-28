
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of an harmonic oscillator using ODE's methods.
// last updated:    28/08/2022


#include "physim.hpp"
#include "gplot++.h"


using namespace physics; 
using namespace measurements; 
using namespace objects; 
using namespace oscillators; 


int main() {

    Gnuplot plot{};
    std::vector<double> t{}, x{};
    unsigned int count{};
    measurement h(1.e-4, s);  
    objects::time tmax(70, s); 
    harmonic oscillator(measurement(1, s.inv())); 
    mass mas(measurement(3, kg), position(0 * m, 0 * m, 0 * m), velocity(1 * mps, 0 * mps, 0 * mps));
    timer timer; 

    plot.redirect_to_png("images/harmonic.png");
    plot.set_xlabel("Time [s]");
    plot.set_ylabel("Position [m]"); 
    
    std::cout << "\nSimulation of an harmonic oscillator with Euler's method\n";
    std::cout << "duration = ";
    tmax.print_time(); 
    std::cout << "increment = "; 
    h.print_measurement(); 
    std::cout << "\nInitial conditions:\n"; 
    oscillator.print_omega(); 
    mas.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 

    // Euler's method
    timer.start(); 
    while (oscillator.get_time() < tmax.get_time()) {
        if (count % 100 == 0) {
            t.emplace_back(oscillator.get_time().value()); 
            x.emplace_back(mas.get_x().value());
        }
        mas.set_pos_vel(oscillator.euler(mas.get_pos_vel(), h)); 
        count++; 
    }
    timer.pause(); 

    std::cout << "\nSimulation with Euler ended\n"; 
    timer.print();    
    std::cout << "\nFinal conditions:\n"; 
    mas.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 
    std::cout << "\n"; 

    plot.plot(t, x, "Euler"); 
    std::cout << "\nPlotting ended (1/2) \n"; 

    mas.set_pos_vel(position(0., 0., 0.), velocity(1., 0., 0.)); 
    oscillator.reset_time(); 
    t.clear(); 
    x.clear(); 
    count = 0;

    std::cout << "\nSimulation of an harmonic oscillator with RK4's method for ";
    tmax.print_time(); 
    std::cout << "\nInitial conditions:\n"; 
    oscillator.print_omega(); 
    mas.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 
    std::cout << "\n"; 

    // Runge Kutta's method
    while (oscillator.get_time() < tmax.get_time()) {
        if (count % 100 == 0) {
            t.emplace_back(oscillator.get_time().value()); 
            x.emplace_back(mas.get_x().value());
        }
        mas.set_pos_vel(oscillator.euler(mas.get_pos_vel(), h)); 
        count++; 
    }
    timer.pause(); 

    std::cout << "\nSimulation with RK4 ended\n"; 
    timer.print();    
    std::cout << "\nFinal conditions:\n"; 
    mas.print_pos_vel(); 
    oscillator.print_time(); 
    std::cout << "\n"; 

    plot.plot(t, x, "RK4"); 
    plot.show(); 
    std::cout << "\nPlotting ended (2/2) \n"; 

    return 0; 

}
