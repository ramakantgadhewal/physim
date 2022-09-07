
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of a forced_damped oscillator using ODE's methods.
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
    objects::time tmax(100, s); 
    forced_damped oscillator(measurement(10, s.inv()), measurement(15, s.inv()), measurement(1. / 30., s.inv())); 
    mass mas(measurement(3, kg), position(0 * m, 0 * m, 0 * m), velocity(0 * mps, 0 * mps, 0 * mps));
    timer timer; 

    plot.redirect_to_png("images/forced_damped.png");
    plot.set_xlabel("Time [s]");
    plot.set_ylabel("Position [m]"); 
    
    std::cout << "\nSimulation of a forced-damped oscillator with RK4's method\n";
    std::cout << "duration = ";
    tmax.print_time(); 
    std::cout << "increment = "; 
    h.print_measurement(); 
    std::cout << "\nInitial conditions:\n"; 
    oscillator.print_omega0(); 
    oscillator.print_omega1(); 
    oscillator.print_alpha(); 
    mas.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 

    // RK4's method
    timer.start(); 
    while (oscillator.get_time() < tmax.get_time()) {
        if (count % 100 == 0) {
            t.emplace_back(oscillator.get_time().value()); 
            x.emplace_back(mas.get_x().value());
        }
        mas.set_pos_vel(oscillator.rk4(mas.get_pos_vel(), h)); 
        count++; 
    }
    timer.pause(); 

    std::cout << "\nSimulation with RK4 ended\n"; 
    timer.print();    
    std::cout << "\nFinal conditions:\n"; 
    mas.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 
    std::cout << "\n"; 

    plot.plot(t, x, "RK4"); 
    plot.show();
    std::cout << "\nPlotting ended (1/1) \n"; 

    return 0; 

}
