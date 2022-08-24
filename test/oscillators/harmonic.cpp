
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of an harmonic oscillator using ODE's methods.
// last updated:    24/08/2022


#include "physim.h"
#include "gplot++.h"


using namespace physim; 
using namespace physics; 
using namespace objects; 


int main() {

    Gnuplot plot{};
    plot.redirect_to_png("images/harmonic.png");
    plot.set_xlabel("Time [s]");
    plot.set_ylabel("Position [m]"); 

    std::vector<double> t, x;
    unsigned int count{};
    
    oscillators::harmonic oscillator(measurements::fixed_measurement(1, units::rad / units::s)); 
    mass m(measurements::fixed_measurement(3, units::kg), position(0., 0., 0.), velocity(1., 0., 0.));
    objects::time tmax(70, units::s); 
    timer timer; 
        
    std::cout << "\nSimulation of an harmonic oscillator with Euler's method for ";
    tmax.print_time(); 
    std::cout << "\nInitial conditions:\n"; 
    oscillator.print_omega(); 
    m.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 
    std::cout << "\n"; 

    // Euler's method
    timer.start(); 
    while (oscillator.get_time().value() < tmax.get_time().value()) {
        if (count % 100 == 0) {
            t.emplace_back(oscillator.get_time().value()); 
            x.emplace_back(m.get_x().value());
        }
        m.set_pos_vel(oscillator.euler(m.get_pos_vel(), 1.e-5)); 
        count++; 
    }
    timer.pause(); 

    std::cout << "\nSimulation with Euler ended\n"; 
    timer.print();    
    std::cout << "\nFinal conditions:\n"; 
    m.print_pos_vel(); 
    oscillator.print_time(); 
    std::cout << "\n"; 

    plot.plot(t, x, "Euler"); 
    std::cout << "\nPlotting ended (1/2) \n"; 

    m.set_pos_vel(position(0., 0., 0.), velocity(1., 0., 0.)); 
    oscillator.reset_time(); 
    t.clear(); 
    x.clear(); 
    count = 0;

    std::cout << "\nSimulation of an harmonic oscillator with RK4's method for ";
    tmax.print_time(); 
    std::cout << "\nInitial conditions:\n"; 
    oscillator.print_omega(); 
    m.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 
    std::cout << "\n"; 

    // Runge Kutta's method
    timer.start(); 
    while (oscillator.get_time().value() < tmax.get_time().value()) {
        if (count % 100 == 0) {
            t.emplace_back(oscillator.get_time().value()); 
            x.emplace_back(m.get_x().value());
        }
        m.set_pos_vel(oscillator.euler(m.get_pos_vel(), 1.e-5)); 
        count++; 
    }

    std::cout << "\nSimulation with RK4 ended\n"; 
    timer.print();    
    std::cout << "\nFinal conditions:\n"; 
    m.print_pos_vel(); 
    oscillator.print_time(); 
    std::cout << "\n"; 

    plot.plot(t, x, "RK4"); 
    plot.show(); 
    std::cout << "\nPlotting ended (2/2) \n"; 

    return 0; 

}
