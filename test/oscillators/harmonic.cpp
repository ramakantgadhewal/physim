
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
    
    oscillators::harmonic oscillator(measurements::fixed_measurement(10, units::rad / units::s)); 
    mass m(measurements::fixed_measurement(3, units::kg), position(3., 3., 3.), velocity(1., 0., 0.));
    m.set_position(m.get_position() / 3); 
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
    // while (oscillator.get_time() < tmax.get_time()) {
        t.emplace_back(oscillator.get_time().value()); 
        x.emplace_back(m.get_x().value());
        m.set_pos_vel(oscillator.euler(m.get_pos_vel())); 

    // }

    // m.print_pos_vel(); 
    timer.pause(); 
    std::cout << "\nSimulation ended\n"; 
    timer.print();    
    std::cout << "Final conditions:\n"; 
    m.print_pos_vel(); 

    plot.plot(t, x, "Euler"); 
    plot.show(); 


    // pos.set_position({0., 0., 0.}, {1., 0., 0.}); 
    // oscill.reset_time(); 
    // time.clear(); 
    // coord_x.clear(); 
    
    // // Runge Kutta's method
    // start = std::chrono::system_clock::now();
    // while (oscill.get_time() < tmax) {
    //     time.push_back(oscill.get_time()); 
    //     coord_x.push_back(pos.get_coord_x());
    //     pos.set_position(oscill.rk4(pos.get_position(), h)); 
    //     oscill.increase_time(h);
    // }    
    // end = std::chrono::system_clock::now();
    // elapsed_seconds = end - start;
    // end_time = std::chrono::system_clock::to_time_t(end);
    
    // std::cout << "Simulation with Runge Kutta's method ended\n"; 
    // std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n" << std::endl;
    
    // plot.plot(time, coord_x, "RK4"); 
    // plot.show(); 
    // std::cout << "Plot ended" << std::endl; 

    return 0; 

}
