

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral class
// last updated:    28/08/2022


#include "physim.hpp"
#include "gplot++.h"


using namespace math; 
using namespace tools; 
using namespace functions; 
using namespace physics::objects;


void test_midpoint(integral& integral, const sine& sin) {

    integral.midpoint(0., constants::pi, sin, 10);
    assert(are_close(integral.value(), 2.0082484079079745, 1.e-16));
            
    integral.midpoint(0., constants::pi, sin, 100);
    assert(are_close(integral.value(), 2.000082249070986, 1.e-15));
            
    integral.midpoint(constants::pi, 0., sin, 10);
    assert(are_close(integral.value(), -2.0082484079079745, 1.e-17));
    
    integral.midpoint(0., 1., sin, 10);
    assert(are_close(integral.value(), 0.45988929071851814, 1.e-17));
    
    integral.midpoint(1., 2., sin, 30);
    assert(are_close(integral.value(), 0.9564934239032155, 1.e-15));
     
    std::cout << "midpoint test passed successfully \n"; 

}


void test_trapexoid(integral& integral, const sine& sin) {

    integral.simpson(0., constants::pi, sin, 10); 
    assert(are_close(integral.value(), 2.0001095173150043, 1.e-16));  

    integral.simpson(0., math::constants::pi, sin, 100); 
    assert(are_close(integral.value(), 2.000000010824504, 1.e-15));  

    integral.simpson(0., 1., sin, 10); 
    assert(are_close(integral.value(), 0.45969794982382056, 1.e-17));  

    integral.simpson(1., 2., sin, 30); 
    assert(are_close(integral.value(), 0.9564491489761575, 1.e-16));  
    
    std::cout << "simpson test passed successfully \n"; 

}


void test_simpson(integral& integral, const sine& sin) {

    integral.trapexoid(0., math::constants::pi, sin, 10); 
    assert(are_close(integral.value(), 1.9835235375094546, 1.e-16));  

    integral.trapexoid(0., math::constants::pi, sin, 100); 
    assert(are_close(integral.value(), 1.9998355038874436, 1.e-16));  

    integral.trapexoid(0., 1., sin, 10); 
    assert(are_close(integral.value(), 0.45931454885797635, 1.e-17));  

    integral.trapexoid(1., 2., sin, 30); 
    assert(are_close(integral.value(), 0.956360580669458, 1.e-15));  

    std::cout << "trapexoid test passed successfully \n"; 

}


void plot_error_runtime(const std::vector<double>& precisions, integral& integral, const sine& sin, const char* title) {

    Gnuplot plot; 
    plot.redirect_to_png(title);
    plot.set_logscale(Gnuplot::AxisScale::LOGXY);
    plot.set_xlabel("precision");
    plot.set_ylabel("error"); 
    std::vector<double> error_midpoint{}, error_trapexoid{}, error_simpson{}; 
    timer timer;

    timer.start();
    for (const auto& precision : precisions) {
        std::cout << "precision: " << precision << "\n";
        integral.midpoint_fixed(0., constants::pi, sin, precision);
        error_midpoint.emplace_back(integral.error());
        integral.trapexoid_fixed(0., constants::pi, sin, precision);
        error_trapexoid.emplace_back(integral.error());        
        integral.simpson_fixed(0., constants::pi, sin, precision);
        error_simpson.emplace_back(integral.error());   
    }
    timer.pause();
    std::cout << "\ncalculations are finished \n";
    timer.print();

    timer.start();
    plot.plot(precisions, precisions, "worst case", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_midpoint, "midpoint", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_trapexoid, "trapexoid", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_simpson, "simpson", Gnuplot::LineStyle::LINESPOINTS);
    plot.show();
    timer.pause();
    std::cout << "\nplotting is finished \n";
    timer.print();

}


void plot_error_known_value(const std::vector<double>& precisions, integral& integral, const sine& sin, const char* title) {

    Gnuplot plot; 
    plot.redirect_to_png(title);
    plot.set_logscale(Gnuplot::AxisScale::LOGXY);
    plot.set_xlabel("precision");
    plot.set_ylabel("error"); 
    std::vector<double> error_midpoint{}, error_trapexoid{}, error_simpson{}; 
    timer timer;

    timer.start();
    for (const auto& precision : precisions) {
        std::cout << "precision: " << precision << "\n";
        integral.midpoint_fixed(0., constants::pi, sin, precision);
        error_midpoint.emplace_back(std::fabs(integral.value() - 2.0));
        integral.trapexoid_fixed(0., constants::pi, sin, precision);
        error_trapexoid.emplace_back(std::fabs(integral.value() - 2.0));        
        integral.simpson_fixed(0., constants::pi, sin, precision);
        error_simpson.emplace_back(std::fabs(integral.value() - 2.0));   
    }
    timer.pause();
    std::cout << "\ncalculations are finished \n";
    timer.print();

    timer.start();
    plot.plot(precisions, precisions, "worst case", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_midpoint, "midpoint", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_trapexoid, "trapexoid", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_simpson, "simpson", Gnuplot::LineStyle::LINESPOINTS);
    plot.show();
    timer.pause();
    std::cout << "\nplotting is finished \n";
    timer.print();

}


int main() {

    integral integral; 
    sine sin;
    timer timer; 

    test_midpoint(integral, sin);
    test_trapexoid(integral, sin);
    test_simpson(integral, sin);

    std::vector<double> precisions = {1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9, 1.e-10, 1.e-11, 1.e-12};

    plot_error_known_value(precisions, integral, sin, "images/error_known_value.png");
    plot_error_runtime(precisions, integral, sin, "images/error_runtime.png");

    return 0;

}