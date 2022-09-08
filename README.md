# physim
**physim** is a c++ _header-only_ library for computational physics.
The goal of this project is to offer the user the opportunity to do his physics' things without having too much trouble with the c++ language and syntax. 


# Table of Contents
* [How to install](#how-to-install)
* [namespace math](#namespace-math)
    * [descriptive_statistics](#descriptive_statistics)
    * [functions](#functions)
    * [tools](#tools)
        * [integral](#integrals)
        * [random_generator](#random_generator)
        * [ode_solver](#ode_solver)
* [namespace physics](#namespace-physics)
    * [constants](#constants)
    * [units](#units)
    * [measurements](#measurements)
    * [position and velocity](#position-and-velocity)
    * [objects](#objects)
        * [material point](#material-point)
        * [mass and charge](#mass-and-charge)
        * [particle](#particle)

# How to install
Download the header file [physim.hpp](https://github.com/lorenzoliuzzo/physim/blob/e0432f73e1ba4ade984c00e8e4b08537f8b42e27/physim.hpp) from terminal typing the following command 
```
wget https://raw.githubusercontent.com/lorenzoliuzzo/physim/master/physim.hpp
``` 

Once you have downloaded the file, you can include it in your desireded .cpp file as 
``` c++
#include "physim.hpp"
```

# namespace math
In this namespace there are defined a few basic tools that allows you to do your calculations. 

## descriptive_statistics
``` c++
using namespace math::descriptive_statistics; 

int main() {

    std::vector<double> data = utilities::read_from_file<double>("data.dat"); 
    double mean = mean(data); 
    double variance = variance(data); 
    double standard_deviation = sd(data); 
    double sd_mean = sdom(data); 
    double median = median(data);
    
    double expected_value = 3.45; 
    double gdl = 3; 
    double chi_squared = chi_sq(data, expected_value);
    double chi_squared_reduced = chi_sq_r(data, expected_value, gdl);
    
    return 0; 
}
```

## functions
In this namespace there are defined some of the most common functions that can be used alone or combined toghether using a functor, providing the two functions and the type of operation (```'+', '-', '*', '/', '^', 'c'```). 
``` c++
using namespace math::functions; 
using namespace math::constants; 

int main() {

    sine sine(3, 2); // y = 3 * sin(2 * x)
    logarithm log(e, 2, 3); // y = 2 * log_5 (3 * x)
    cube cb(4, -2, 0, 1); // y = 4 * x^3 - 2 * x^2 + 1
    square_root sq_rt; // y = sqrt(x)

    // you can print the equation
    sine.print_equation(); 

    // you can compone functions using functors
    functor f('+', sine, cb); // y = 3 * sin(2 * x) + 4 * x^3 - 2 * x^2 + 1
    functor g('c', sq_rt, sine); // y = sqrt(3 * sin(2 * x))

    // you can use funtors as functions as well
    functor h('/', f, g); // y = f / g

    // you can evaluate the function in a specific point and print the value
    double value = h.eval(pi); 
    std::cout << value << "\n"; // old school
    h.print_eval(pi); // cool way

    return 0; 
}
```
The available functions are ```line, quadratic, cubic, square_root, cubic_root, exponential, logarithm, sine, cosine, tangent```. 
By default, the multiplicative parameters are fixed as 1, the adding parameters are fixed as 0 and the basis parameters for the exponential and the logarithm are fixed as ```math::constants::e```. 


## tools

### integrals
``` c++
using namespace math::tools; 
using namespace math::constants; 

int main() {

    integral integral; 
    functions::sine sine; 
    unsigned int steps = 100; 
    double precision = 1.e-6;
    double max_value = 1; 

    // integrate with the midpoint method (default steps = 100) 
    integral.midpoint(0, pi, sine, steps); 

    double value = integral.value(); // get the integral value
    integral.print_value(1.e-6); // print the value

    // integrate with a fixed level of error (default precision = 1.e-6)
    integral.midpoint_fixed(0, pi, sine, precision);

    double error = integral.error(); // get the integral error
    integral.print_error(); // print the error (only after a fixed precision method)
    integral.print_integral(); // print the value and the error


    return 0; 
}
```
The available integration methods are ```midpoint, trapexoid, simpson, mean,  hit_or_miss```, in both "fixed steps" and "fixed precision" (```_fixed```) mode.


### random_generator
``` c++
using namespace math::tools; 

int main() {

    random_generator rg; 
    double mean = 2; 
    double sigma = 0.1; 

    double rand = rg.rand(3, 6); // default are 0 and 1
    double exp = rg.exp(mean); 
    double gauss_BM = rg.gauss_box_muller(mean, sigma); 
    double gauss_AR = rg.gauss_accept_reject(mean, sigma); 

    rg.seed(34); // setting the seed
    rg.seed(); // getting the seed

    return 0; 
}
```

### ode_solver
In this virtual class are defined some of the method for evaluating ordinary differential equations. The available methods are ```euler, rk4```.


# namespace physics

## constants

## units

## measurements

## position & velocity

## objects

