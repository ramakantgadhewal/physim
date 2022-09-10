

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::integral errors
// last updated:    10/09/2022


#include "physim.hpp"
#include "gplot++.h"


using namespace math::tools; 
using namespace math::functions; 
using namespace math::constants; 

using physics::objects::timer;


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
        integral.midpoint_fixed(0., pi, sin, precision);
        error_midpoint.emplace_back(integral.error());
        integral.trapexoid_fixed(0., pi, sin, precision);
        error_trapexoid.emplace_back(integral.error());        
        integral.simpson_fixed(0., pi, sin, precision); 
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
        integral.midpoint_fixed(0., pi, sin, precision);
        error_midpoint.emplace_back(std::fabs(integral.value() - 2.0));
        integral.trapexoid_fixed(0., pi, sin, precision);
        error_trapexoid.emplace_back(std::fabs(integral.value() - 2.0));        
        integral.simpson_fixed(0., pi, sin, precision);
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


void plot_error(const std::vector<double>& precisions, integral& integral, const sine& sin, const char* title) {
    
    Gnuplot plot; 
    plot.redirect_to_png(title);
    plot.set_logscale(Gnuplot::AxisScale::LOGXY);
    plot.set_xlabel("precision");
    plot.set_ylabel("error"); 
    std::vector<double> error_runtime_midpoint{}, error_runtime_trapexoid{}, error_runtime_simpson{}; 
    std::vector<double> error_known_value_midpoint{}, error_known_value_trapexoid{},  error_known_value_simpson{}; 
    timer timer;

    timer.start();
    for (const auto& precision : precisions) {
        integral.midpoint_fixed(0., pi, sin, precision);
        error_known_value_midpoint.emplace_back(std::fabs(integral.value() - 2.0));
        error_runtime_midpoint.emplace_back(integral.error());
        integral.trapexoid_fixed(0., pi, sin, precision);
        error_known_value_trapexoid.emplace_back(std::fabs(integral.value() - 2.0));
        error_runtime_trapexoid.emplace_back(integral.error());        
        integral.simpson_fixed(0., pi, sin, precision);
        error_known_value_simpson.emplace_back(std::fabs(integral.value() - 2.0));
        error_runtime_simpson.emplace_back(integral.error());
    }
    timer.pause();
    std::cout << "\ncalculations are finished \n";
    timer.print();

    timer.start();
    plot.plot(precisions, precisions, "worst case", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_runtime_midpoint, "rt midpoint", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_known_value_midpoint, "kv midpoint", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_runtime_trapexoid, "rt trapexoid", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_known_value_trapexoid, "kv trapexoid", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_runtime_simpson, "rt simpson", Gnuplot::LineStyle::LINESPOINTS);
    plot.plot(precisions, error_known_value_simpson, "kv simpson", Gnuplot::LineStyle::LINESPOINTS);
    plot.show(); 
    timer.pause();
    std::cout << "\nplotting is finished \n";
    timer.print();

}


int main() {

    integral integral; 
    sine sin;
    std::vector<double> precisions = {1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9, 1.e-10, 1.e-11, 1.e-12};

    plot_error_known_value(precisions, integral, sin, "images/error_known_value.png");
    plot_error_runtime(precisions, integral, sin, "images/error_runtime.png");
    plot_error(precisions, integral, sin, "images/error.png");

    return 0;

} 

