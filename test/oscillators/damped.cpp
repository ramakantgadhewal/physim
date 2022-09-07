
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     Plot of a damped oscillator using ODE's methods.
// last updated:    28/08/2022


#include "physim.hpp"
#include "gplot++.h"


using namespace physics; 
using namespace measurements; 
using namespace objects; 
using namespace oscillators; 


void plot_damped_oscillator(const measurement& omega, 
                            const measurement& alpha, 
                            const material_point& mp,
                            const objects::time& time_max, 
                            const measurement& increment,
                            const char* filename_and_type) {

    std::vector<double> t{}, x{};
    // unsigned int count{};
    timer timer; 
    damped oscillator(omega, alpha); 
    material_point mass(mp.get_pos_vel());
    Gnuplot plot{};
    plot.redirect_to_png(filename_and_type);
    plot.set_xlabel("Time [s]");
    plot.set_ylabel("Position [m]"); 
    
    std::cout << "\nSimulation of a damped oscillator with RK4's method\n";
    std::cout << "duration = ";
    time_max.print_time(); 
    std::cout << "increment = "; 
    increment.print_measurement(); 
    std::cout << "\nInitial conditions:\n"; 
    oscillator.print_omega(); 
    oscillator.print_alpha(); 
    mass.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 

    // RK4's method
    timer.start(); 
    while (oscillator.get_time() < time_max.get_time()) {
        // if (count % 100 == 0) {
            t.emplace_back(oscillator.get_time().value()); 
            x.emplace_back(mass.get_x().value());
        // }
        mass.set_pos_vel(oscillator.rk4(mass.get_pos_vel(), increment)); 
        // count++; 
    }
    timer.pause(); 

    std::cout << "\nSimulation with RK4 ended\n"; 
    timer.print();    
    std::cout << "\nFinal conditions:\n"; 
    mass.print_pos_vel(); 
    std::cout << "time = ";
    oscillator.print_time(); 
    std::cout << "\n"; 

    plot.plot(t, x, "RK4"); 
    plot.show();
    std::cout << "\nPlotting ended (1/1) \n"; 

}

int main() {

    measurement omega(40, s.inv()), alpha(3, s.inv());
    material_point mass(position(3 * m, 0 * m, 0 * m), velocity(1 * mps, 0 * mps, 0 * mps));
    objects::time time_max(50, s); 
    measurement increment(1.e-4, s); 

    // plotting the damped oscillator
    plot_damped_oscillator(omega, alpha, mass, time_max, increment, "images/damped.png");


    return 0; 

}
