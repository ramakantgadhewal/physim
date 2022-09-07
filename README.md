# physim
**physim** is a c++ header-only namespace for computational physics.
The goal of this project is to offer the user the opportunity to do his physics' things without having too much trouble with the c++ language and syntax. 


# Table of Contents
* [How to install and include it](#how_to_install_and_include_it)
* [namespace math](#namespace-math)
  * [constants](#constants)
  * [algebra](#algebra)
  * [descriptive_statistics](#descriptive_statistics)
  * [functions](#functions)
  * [integrals](#integrals)
  * [random_generator](#random_generator)
  * [ode_solver](#ode_solver)
* [namespace physics](#namespace-physics)


# How to install and include it
Download the header file [physim.hpp](https://github.com/lorenzoliuzzo/physim/blob/e0432f73e1ba4ade984c00e8e4b08537f8b42e27/physim.hpp) and include it in your .cpp as 
``` c++
#include "physim.hpp"
```


# namespace math
In this namespace there are defined a few basic tools that allows you to do your calculations. 

## constants
Here there are defined some useful math constants that can be used as constexpr as shown by the example below. 
``` c++
using namespace math::constants; 

int main() {

   double r = 1; 
   double cfr = 2 * pi * r; 
   
   return 0;
}
```

## algebra

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

int main() {
   sine sin(3, 2); // y = 3 * sin(2 * x)
   cube cb(4, -2, 0, 1); // y = 4 * x^3 - 2 * x^2 + 1
   square_root sqrt; // y = sqrt(x)
   
   // you can print the equation
   sin.print_equation(); 
   
   // you can compone functions using functors
   functor f('+', sin, cb); // y = 3 * sin(2 * x) + 4 * x^3 - 2 * x^2 + 1
   functor g('c', sqrt, sin); // y = sqrt(3 * sin(2 * x))
   
   // you can use funtors as functions as well
   functor h('/', f, g); // y = f / g
   
   // you can evaluate the function in a specific point and print the value
   double value = h.eval(constants::pi); 
   std::cout << value << "\n"; // old school
   h.print_eval(constants::pi); // cool way
   
   return 0; 
}
```

## integrals

## random_generator

## ode_solver

# namespace physics

## constants

## units

## measurements

## position & velocity

## objects

