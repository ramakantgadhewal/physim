
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     testing the math::tools::random_generator class
// last updated:    20/08/2022

#include "physim.hpp"
#include "gplot++.h"


using namespace math::tools; 


void plot_unif(const char* title, random_generator& rg, unsigned int& nmax, const double& a = 0., const double& b = 1., const unsigned int& bins = 30) {

    Gnuplot plot{};
    plot.redirect_to_png(title);
    plot.set_xlabel("Values");
    plot.set_ylabel("Number of counts");
    plot.set_yrange(0., NAN); 

    std::vector<double> rand; 
    for (unsigned int i{}; i < nmax; i++) rand.emplace_back(rg.unif(a, b));
    plot.histogram(rand, bins, "rand");
    plot.show();
    std::cout << "plotting ended\n"; 

}

void plot_exp(const char* title, random_generator& rg, unsigned int& nmax, const double& mean, const unsigned int& bins = 10) {

    Gnuplot plot{};
    plot.redirect_to_png(title);
    plot.set_xlabel("Values");
    plot.set_ylabel("Number of counts");
    plot.set_yrange(0., NAN); 

    std::vector<double> exp; 
    for (unsigned int i{}; i < nmax; i++) exp.emplace_back(rg.exp(mean));
    plot.histogram(exp, bins, "exponential");
    plot.show();
    std::cout << "plotting ended\n"; 

}

void plot_gauss_bm(const char* title, random_generator& rg, unsigned int& nmax, const double& mean, const double& sigma, const unsigned int& bins = 10) {

    Gnuplot plot{};
    plot.redirect_to_png(title);
    plot.set_xlabel("Values");
    plot.set_ylabel("Number of counts");
    plot.set_yrange(0., NAN); 

    std::vector<double> gauss; 
    for (unsigned int i{}; i < nmax; i++) 
        gauss.emplace_back(rg.gauss_box_muller(mean, sigma));
    plot.histogram(gauss, bins, "gauss Box-Muller");
    plot.show();
    std::cout << "plotting ended\n"; 
    
}

void plot_gauss_ar(const char* title, random_generator& rg, unsigned int& nmax, const double& mean, const double& sigma, const unsigned int& bins = 10) {

    Gnuplot plot{};
    plot.redirect_to_png(title);
    plot.set_xlabel("Values");
    plot.set_ylabel("Number of counts");
    plot.set_yrange(0., NAN); 

    std::vector<double> gauss; 
    for (unsigned int i{}; i < nmax; i++) 
        gauss.emplace_back(rg.gauss_accept_reject(mean, sigma));
    plot.histogram(gauss, bins, "gauss accept-reject");
    plot.show();
    std::cout << "plotting ended\n"; 
    
}

void plot_lorentzian( const char* title, random_generator& rg, unsigned int& nmax, const double& mean, const double& gamma, const unsigned int& bins = 10) {

    Gnuplot plot{};
    plot.redirect_to_png(title);
    plot.set_xlabel("Values");
    plot.set_ylabel("Number of counts");
    plot.set_yrange(0., NAN); 
    // plot.set_xrange(-10, 10); 

    std::vector<double> lor; 
    for (unsigned int i{}; i < nmax; i++) {
        double x = rg.lorentzian(mean, gamma); 
        while (x < mean - 10 * gamma || x > mean + 10 * gamma) x = rg.lorentzian(mean, gamma);
        lor.emplace_back(x); 
    }
    std::cout << *(std::max_element(lor.begin(), lor.end())) << "\n";
    std::cout << *(std::min_element(lor.begin(), lor.end())) << "\n";

    // std::cout << "min = " << min_iter << "\n"; 

    plot.histogram(lor, bins, "lorentzian");
    plot.show();
    std::cout << "plotting ended\n"; 
    
}


int main() {    

    random_generator rg;
    rg.set_up(); 
    unsigned int nmax{10000000}; 
    // plot_unif("images/unif.png", rg, nmax, 3., 6., 50);
    // plot_exp("images/exp.png", rg, nmax, 3, 50); 
    // plot_gauss_bm("images/gauss_bm.png", rg, nmax, 2., 0.2, 50); 
    // plot_gauss_ar("images/gauss_ar.png", rg, nmax, 2., 0.2, 50); 
    plot_lorentzian("images/lorentzian.png", rg, nmax, 2., 0.2, 50);  
    
    return 0; 
}