

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physim(namespace) containing the basic tools for computational physics. 
// last updated:    10/09/2022


#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <functional>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>


namespace utilities {

    // reading the template data from a file to an std::vector<>
    template <typename T>
    std::vector<T> read_from_file(const char* filename) {
        std::vector<T> vec; 
        std::ifstream filein; 
        filein.open(filename); 
        T appo; 
        if (filein.fail()) {
            std::cerr << "Error! It was not possible to open the selected file \n"; 
            exit(11); 
        } 
        for (unsigned int i{}; !filein.eof(); i++) {
            filein >> appo;
            vec.emplace_back(appo); 
        }
        filein.close();
        return vec;    
    }

} // namespace utilities


namespace math {
            

    // namespace defining some usefull constants
    namespace constants {

        constexpr int32_t max_neg(uint32_t n_bits) { return -(int32_t(1U << (n_bits - 1))); }
        constexpr double infinity = std::numeric_limits<double>::infinity();
        constexpr double invalid_conversion = std::numeric_limits<double>::signaling_NaN();

        constexpr double pi = 3.14159265358979323846;
        constexpr double e = 2.7182818284590452353603;

    } // namespace constants
    

    // namespace defining some usefull operation
    namespace algebra {

        // generate the square power of a value
        template<typename T> 
        constexpr T square(const T& value) { return std::pow(value, 2); }

        // generate the cubic power of a value
        template<typename T> 
        constexpr T cube(const T& value) { return std::pow(value, 3); }

        // generate small integer powers of a value (1, 0, -1)
        template<typename T> 
        constexpr T pow_small(const T& value, const int& power) { 
            return (power == 1) ? value : ((power == -1) ? T(1.0) / value : T(1.0));
        }

        // generate root power of a value
        constexpr double root(const double& value, const int& power) {
            switch (power) {
                case 0: 
                    return 1.0;
                case 1: 
                    return value;
                case -1: 
                    return 1.0 / value;
                case 2: 
                    if (value < 0.0) { return constants::invalid_conversion; }
                    else return std::sqrt(value);
                case -2:
                    if (value < 0.0) { return constants::invalid_conversion; }
                    else return std::sqrt(1.0 / value);
                case 3: 
                    return std::cbrt(value);
                case -3: 
                    return std::cbrt(1.0 / value);
                case 4: 
                    if (value < 0.0) { return constants::invalid_conversion; }
                    else return std::sqrt(std::sqrt(value));
                case -4: 
                    if (value < 0.0) { return constants::invalid_conversion; }
                    else return std::sqrt(std::sqrt(1.0 / value));
                default:
                    if (value < 0.0 && power % 2 == 0) { return constants::invalid_conversion; }
                    else return std::pow(value, 1.0 / power);
            }
        }
        
        // generate the square root power of a value
        constexpr double sqrt(const double& value) { return root(value, 2); }

        // generate the cubic root power of a value
        constexpr double cbrt(const double& value) { return root(value, 3); }


    } // namespace algebra


    // namespace defining some descriptive statistic functions for std::vector<>
    namespace descriptive_statistics {

        // mean of an std::vector
        template <typename K>
        double mean(const std::vector<K>& v) {
            if (v.size() == 0) return 0.0;
            double accu{};
            for (size_t i{}; i < v.size(); i++) {
                accu = static_cast<double>(i) / static_cast<double>(i+1) * accu + 1.0 / static_cast<double>(i+1) * v[i];
            }
            return accu;
        }

        // median of an std::vector
        template <typename K>
        double median(std::vector<K>& v) {
            if (std::is_sorted(v.begin(), v.end()) == false) std::sort(v.begin(), v.end());
            if (v.size() % 2 != 0) return v[v.size() / 2];
            else return (v[v.size() / 2] + v[(v.size() / 2) - 1]) / 2; 
        }

        // variance of an std::vector
        template <typename K>
        double variance(const std::vector<K>& v) {
            if(v.size() == 0) return 0.0;
            double result{};
            double old_average{};
            double average{};
            for(size_t i{}; i < v.size(); i++){
                old_average = average;
                average = static_cast<double>(i) / static_cast<double>(i+1) * average + 1.0 / static_cast<double>(i+1) * v[i];
                result = 1.0 / static_cast<double>(i+1) * (static_cast<double>(i) * result + algebra::square(v[i]) + static_cast<double>(i) * algebra::square(old_average)) - algebra::square(average);
            }
            return result;
        }

        // standard deviation of an std::vector
        template <typename K>
        double sd(const std::vector<K>& v) {
            return algebra::sqrt(variance(v));
        }

        // standard deviation of mean of an std::vector
        template <typename K>
        double sdom(const std::vector<K>& v) {
            return sd(v) / algebra::sqrt(v.size());
        }

        // chi squared of an std::vector
        template <typename K>
        double chi_sq(const std::vector<K>& v, const double& expected_value) {
            double accu{}; 
            for (auto x : v) accu += std::pow(x - expected_value, 2) / sd(v); 
            return accu; 
        }           

        // chi squared reduced of an std::vector
        template <typename K>
        double chi_sq_r(const std::vector<K>& v, const double& expected_value, const unsigned int& gdl) {
            double accu{}; 
            for (auto x : v) accu += std::pow(x - expected_value, 2) / sd(v); 
            return accu / gdl; 
        }
        

    } // namespace descriptive_statistics


    // namespace defining some of the most common math functions
    namespace functions {

        // function_base 
        class function_base {

            public: 
            
                // =============================================
                // virtual destructor
                // =============================================     
            
                virtual ~function_base() {}

            
                // =============================================
                // eval methods
                // =============================================
            
                virtual constexpr double eval(const double& x) const = 0; 
            
                constexpr int signum(const double& x) const { return (eval(x) == 0. ? 0. : (eval(x) > 0 ? 1. : -1)); }


                // =============================================
                // print methods
                // =============================================
            
                void print_eval(const double& x, const double& precision = 1.e-6) const { std::cout << "f(" << x << ") = " << std::setprecision(precision) << eval(x) << "\n"; }
            
                virtual void print_equation() const = 0;
            

        }; // class function_base
        

        // functor class for composing single functions into some more complex expressions
        class functor : public function_base {

            private:

                // =============================================
                // class members
                // =============================================
            
                char op_; 

                function_base * f_; 

                function_base * g_;
            
            
            public:

                // =============================================
                // constructor and destructor
                // =============================================     
            
                functor(const char& op, function_base* f, function_base* g) noexcept {
                    if (op != '+' && op != '-' && op != '*' && op != '/' && op != '^' && op != 'c') std::cerr << "Wrong functor operation. The possibles operations are ['+', '-', '*', '/', '^', 'c'] \n"; 
                    op_ = op;
                    f_ = f; 
                    g_ = g;
                }
                
                ~functor() {}
            
            
                // =============================================
                // eval methods
                // =============================================
            
                constexpr double eval(const double& x) const override {
                    switch (op_) {
                        case '+': return f_->eval(x) + g_->eval(x);
                        case '-': return f_->eval(x) - g_->eval(x);
                        case '*': return f_->eval(x) * g_->eval(x);
                        case '/': return f_->eval(x) / g_->eval(x);
                        case '^': return std::pow(f_->eval(x), g_->eval(x));
                        case 'c': return f_->eval(g_->eval(x));
                    }
                }

                
                // =============================================
                // print methods
                // =============================================
            
                void print_equation() const override {}
                

        }; // class functor


        // line function
        class line : public function_base {

            private: 

                // =============================================
                // class members
                // =============================================
                
                // y = m * x + q
                double m_, q_; 


            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================     
        
                explicit constexpr line(const double& m = 1, const double& q = 0) noexcept : m_{m}, q_{q} {}
                
                ~line() {}

        
                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void m(const double& m) { m_ = m; } 

                constexpr void q(const double& q) { q_ = q; }
        
                constexpr double m() const { return m_; }

                constexpr double q() const { return q_; } 


                // =============================================
                // eval methods
                // =============================================

                constexpr double eval(const double& x) const override { return m_ * x + q_; }
        
        
                // =============================================
                // print methods
                // =============================================

                void print_equation() const override {
                    std::cout << " y = "; 
                    if (m_ != 1 && m_ != -1) std::cout << m_; 
                    else {
                        if (m_ == 1) std::cout << "+";
                        if (m_ == -1) std::cout << "-";
                    }
                    std::cout << "x"; 
                    if (q_ != 0) {
                        if (q_ > 0) std::cout << " + " << q_ << "\n"; 
                        else std::cout << " - " << std::fabs(q_) << "\n"; 
                    } else std::cout << "\n";
                }  

        
        }; // class line


        // quadratic function
        class quadratic : public function_base {

            private: 

                // =============================================
                // class members
                // =============================================
                
                // y = a * x^2 + b * x + c
                double a_, b_, c_, delta_; 


            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================     
        
                explicit constexpr quadratic(const double& a, const double& b, const double& c) noexcept :
                    a_{a}, b_{b}, c_{c}, delta_{algebra::square(b_) - 4 * a_ * c_} {}
                
                ~quadratic() {}

        
                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void a(const double& a) { a_ = a; } 

                constexpr void b(const double& b) { b_ = b; }

                constexpr void c(const double& c) { c_ = c; }  
        
                constexpr double a() const { return a_; }

                constexpr double b() const { return b_; } 

                constexpr double c() const { return c_; } 

                constexpr double delta() const { return delta_; } 

                // constexpr std::pair<auto, auto> roots() const {
                //     if (delta_ == 0) return std::make_pair(- b_ / (2 * a_), - b_ / (2 * a_));
                //     else if (delta_ > 0) return std::make_pair(- b_ - algebra::sqrt(delta_) / (2 * a_), (- b_ + algebra::sqrt(delta_)) / (2 * a_));
                //     else return std::make_pair<std::complex<double>, std::complex<double>>((- b_, algebra::sqrt(std::fabs(delta_))), (- b_, - algebra::sqrt(std::fabs(delta_))));
                //     return std::make_pair(NAN, NAN); 
                // }


                // =============================================
                // eval methods
                // =============================================

                constexpr double eval_Horner(const double& x) const { return c_ + x * (b_ + x * a_); }

                constexpr double eval(const double& x) const override { return a_ * algebra::square(x) + b_ * x + c_; }
        
        
                // =============================================
                // print methods
                // =============================================

                void print_equation() const override {
                    std::cout << " y = "; 
                    if (a_ != 1 && a_ != -1) std::cout << a_; 
                    else if (a_ == -1) std::cout << "-"; 
                    else std::cout << "x^2 "; 
                    if (b_ != 0) {
                        if (b_ > 0) std::cout << "+ ";
                        else std::cout << "- ";
                        if (b_ != 1 && b_ != -1) std::cout << std::fabs(b_) << "x ";  
                        else std::cout << "x "; 
                    } 
                    if (c_ != 0) {
                        if (c_ > 0) std::cout << "+ ";
                        else std::cout << "- "; 
                        std::cout << std::fabs(c_) << "\n";
                    } else std::cout << "\n";
                }  

        
        }; // class quadratic


        // cubic function
        class cubic : public function_base {
            
            private: 

                // =============================================
                // class members
                // =============================================
                
                // y = a * x^3 + b * x^2 + c * x + d
                double a_, b_, c_, d_; 


            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
            
                explicit constexpr cubic(const double& a, const double& b, const double& c, const double& d) noexcept :
                    a_{a}, b_{b}, c_{c}, d_{d} {}
            
                ~cubic() {}  
            
            
                // =============================================
                // set methods
                // =============================================
        
                constexpr void a(const double& a) { a_ = a; } 

                constexpr void b(const double& b) { b_ = b; }

                constexpr void c(const double& c) { c_ = c; }  
        
                constexpr void d(const double& d) { d_ = d; } 
        
            
                // =============================================
                // get methods
                // =============================================

                constexpr double a() const { return a_; }

                constexpr double b() const { return b_; } 

                constexpr double c() const { return c_; } 

                constexpr double d() const { return d_; } 

            
                // =============================================
                // eval methods
                // =============================================

                constexpr double eval(const double& x) const override { 
                    return a_ * algebra::cube(x) + b_ * algebra::square(x) + c_ * x + d_; 
                }
        
        
                // =============================================
                // print methods
                // =============================================

                void print_equation() const override {
                    std::cout << " y = "; 
                    if (a_ != 1) std::cout << a_; 
                    std::cout << "x^3 "; 
                    if (b_ != 0) {
                        if (b_ != 1) std::cout << "+ " << b_ << "x^2 "; 
                        else std::cout << "+ x^2 "; 
                    }
                    if (c_ != 0) {
                        if (c_ != 1) std::cout << "+ " << c_ << "x"; 
                        else std::cout << "+ x"; 
                    }
                    if (d_ != 0) std::cout << " + " << d_ << "\n";
                    else std::cout << "\n";  
                }

        
        }; // class cubic


        // square_root function
        class square_root : public function_base {
            
            private: 

                // =============================================
                // class members
                // =============================================
                
                // y = c * x ^ (1 / 2)
                double c_; 


            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
            
                explicit constexpr square_root(const double& c = 1) noexcept : c_{c} {}
            
                ~square_root() {}
            
            
                // =============================================
                // set & get methods
                // =============================================
            
                constexpr void c(const double& c) { c_ = c; }

                constexpr double c() const { return c_; }

            
                // =============================================
                // eval methods
                // =============================================

                constexpr double eval(const double& x) const override { return c_ * algebra::sqrt(x); }
                        
        
                // =============================================
                // print methods
                // =============================================

                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c_ != 1) std::cout << c_; 
                    std::cout << "x^(1/2) \n";
                }        
            
            
        }; // class square_root


        // cubic_root function
        class cubic_root : public function_base {
            
            private: 

                // =============================================
                // class members
                // =============================================
                
                // y = c * x ^ (1 / 3)
                double c_; 


            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
            
                explicit constexpr cubic_root(double c = 1) noexcept : c_{c} {}
            
                ~cubic_root() {} 
            
            
                // =============================================
                // set & get methods
                // =============================================
            
                constexpr void c(double c) { c_ = c; }

                constexpr double c() const { return c_; }

            
                // =============================================
                // eval methods
                // =============================================

                constexpr double eval(const double& x) const override { return c_ * algebra::cbrt(x); }
                        
        
                // =============================================
                // print methods
                // =============================================

                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c_ != 1) std::cout << c_; 
                    std::cout << "x^(1/3) \n";
                }        
            

        }; // class cubic_root


        // exponential function
        class exponential : public function_base {
        
            private:
        
                // =============================================
                // class members
                // =============================================
                
                // y = c1 * base ^ (c2 * x)
                double base_, c1_, c2_;

            
            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
        
                explicit constexpr exponential(const double& base = constants::e, const double& c1 = 1, const double& c2 = 1) noexcept :
                    base_{base}, c1_{c1}, c2_{c2} {}

                ~exponential() {} 
        

                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void base(const double& base) { base_ = base; } 
            
                constexpr void c1(const double& c1) { c1_ = c1; }
            
                constexpr void c2(const double& c2) { c2_ = c2; }
        
                constexpr double base() const { return base_; }
        
                constexpr double c1() const { return c1_; }

                constexpr double c2() const { return c2_; }
            
            
                // =============================================
                // eval methods
                // =============================================
        
                constexpr double eval(const double& x) const override { return c1_ * std::pow(base_, c2_ * x); }
                        
            
                // =============================================
                // print methods
                // =============================================
        
                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c1_ != 1) std::cout << c1_;
                    if (base_ != constants::e) std::cout << base_ << "^"; 
                    else std::cout << "e^(";
                    if (c2_ != 1 && c2_ != -1) std::cout << c2_; 
                    if (c2_ == -1) std::cout << "-"; 
                    std::cout << "x) \n"; 
                }


        }; // class exponential


        // logarithm function
        class logarithm : public function_base {
        
            private:
        
                // =============================================
                // class members
                // =============================================
                
                // y = c1 * log(base) (c2 * x) 
                double base_, c1_, c2_; 

            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
        
                explicit constexpr logarithm(const double& base = constants::e, const double& c1 = 1, const double& c2 = 1) noexcept : 
                    base_{base}, c1_{c1}, c2_{c2} {}

                ~logarithm() {} 
        

                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void base(const double& base) { base_ = base; } 
        
                constexpr void c1(const double& c1) { c1_ = c1; }
            
                constexpr void c2(const double& c2) { c2_ = c2; }
                        
                constexpr double base() const { return base_; }

                constexpr double c1() const { return c1_; }

                constexpr double c2() const { return c2_; }
            

                // =============================================
                // eval methods
                // =============================================
        
                constexpr double eval(const double& x) const override { return c1_ * std::log(c2_ * x) / std::log(base_); }
                            
            
                // =============================================
                // print methods
                // =============================================
        
                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c1_ != 1) std::cout << c1_;
                    std::cout << "log";
                    if (base_ != constants::e) std::cout << "_( " << base_ << ")";
                    if (c2_ != 1) std::cout << "(" << c2_ << "x) \n"; 
                    else std::cout << "(x) \n";
                }                


        }; // class logarithm 


        // sine function
        class sine : public function_base {
            
            private: 

                // =============================================
                // class members
                // =============================================   
                
                // y = c1 * sin (c2 * x + c3)
                double c1_, c2_, c3_; 
            
            
            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
        
                explicit constexpr sine(const double& c1 = 1, const double& c2 = 1, const double& c3 = 0.) noexcept : 
                    c1_{c1}, c2_{c2}, c3_{c3} {} 
            
                ~sine() {} 

            
                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void c1(const double& c1) { c1_ = c1; }
            
                constexpr void c2(const double& c2) { c2_ = c2; }

                constexpr void c3(const double& c3) { c3_ = c3; }
            
                constexpr double c1() const { return c1_; }

                constexpr double c2() const { return c2_; }

                constexpr double c3() const { return c3_; }
            
            
                // =============================================
                // eval methods
                // =============================================
        
                constexpr double eval(const double& x) const override { return c1_ * std::sin(c2_ * x); }
                    
            
                // =============================================
                // print methods
                // =============================================
        
                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c1_ != 1) std::cout << c1_;
                    std::cout << "sin(";
                    if (c2_ != 1) std::cout << c2_; 
                    std::cout << "x) \n";
                }


        }; // class sine


        // cosine function
        class cosine : public function_base {
            
            private: 

                // =============================================
                // class members
                // =============================================   
                
                // y = c1 * cos (c2 * x + c3)
                double c1_, c2_, c3_; 
            
            
            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
        
        
                explicit constexpr cosine(const double& c1 = 1, const double& c2 = 1, const double& c3 = 0.) noexcept : 
                    c1_{c1}, c2_{c2}, c3_{c3} {} 
            
                ~cosine() {} 

            
                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void c1(const double& c1) { c1_ = c1; }
            
                constexpr void c2(const double& c2) { c2_ = c2; }

                constexpr void c3(const double& c3) { c3_ = c3; }
            
                constexpr double c1() const { return c1_; }

                constexpr double c2() const { return c2_; }

                constexpr double c3() const { return c3_; }
            

                // =============================================
                // eval methods
                // =============================================
        
                constexpr double eval(const double& x) const override { return c1_ * std::cos(c2_ * x); }
                    
        
                // =============================================
                // print methods
                // =============================================
        
                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c1_ != 1) std::cout << c1_;
                    std::cout << "cos(";
                    if (c2_ != 1) std::cout << c2_; 
                    std::cout << "x) \n";
                }


        }; // class cosine


        // tangent function
        class tangent : public function_base {

            private: 

                // =============================================
                // class members
                // =============================================   
                
                // y = c1 * tan (c2 * x + c3)
                double c1_, c2_, c3_; 
            
            
            public: 
        
                // =============================================
                // constructor and destructor
                // =============================================   
        
                explicit constexpr tangent(const double& c1 = 1, const double& c2 = 1, const double& c3 = 0.) noexcept : 
                    c1_{c1}, c2_{c2}, c3_{c3} {} 
            
                ~tangent() {} 

            
                // =============================================
                // set & get methods
                // =============================================
        
                constexpr void c1(const double& c1) { c1_ = c1; }
            
                constexpr void c2(const double& c2) { c2_ = c2; }

                constexpr void c3(const double& c3) { c3_ = c3; }
            
                constexpr double c1() const { return c1_; }

                constexpr double c2() const { return c2_; }

                constexpr double c3() const { return c3_; }
            

                // =============================================
                // eval methods
                // =============================================
        
                constexpr double eval(const double& x) const override { return c1_ * std::tan(c2_ * x); }
                
        
                // =============================================
                // print methods
                // =============================================

                void print_equation() const override {
                    std::cout << " y = "; 
                    if (c1_ != 1) std::cout << c1_;
                    std::cout << "tan(";
                    if (c2_ != 1) std::cout << c2_; 
                    std::cout << "x) \n";
                }


        }; // class tangent


    } // namespace functions
    

    // namespace defining some usefull operations
    namespace tools {
        

        // check if the two value given differ from each other less than the precision given
        constexpr bool are_close(const double& calculated, const double& expected, const double& epsilon = 1e-6){
            return std::fabs(calculated - expected) <= epsilon;
        }


        // round a value to the expected level of precision of a double
        double cround(double val) {
            std::uint64_t bits;
            std::memcpy(&bits, &val, sizeof(bits));
            bits += 0x800ULL;
            bits &= 0xFFFFFFFFFFFFF000ULL;
            std::memcpy(&val, &bits, sizeof(bits));
            return val;
        }


        // rounding compare for equality on double
        bool compare_round_equals(double val1, double val2) {
            static constexpr double half_precise_precision{5e-13};
            auto v1 = val1 - val2;
            if (v1 == 0.0 || std::fpclassify(v1) == FP_SUBNORMAL) { return true; }
            auto c1 = cround(val1);
            auto c2 = cround(val2);
            return (c1 == c2) ||
                (cround(val2 * (1.0 + half_precise_precision)) == c1) ||
                (cround(val2 * (1.0 - half_precise_precision)) == c1) ||
                (cround(val1 * (1.0 + half_precise_precision)) == c2) ||
                (cround(val1 * (1.0 - half_precise_precision)) == c2);
        }


        // class random_generator for generating pseudo-casual numbers
        class random_generator {

            private: 

                // =============================================
                // class members
                // =============================================        
                
                unsigned int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;


            public: 
        
                // =============================================
                // constructor & destructor
                // =============================================

                random_generator() noexcept {}

                ~random_generator() {}


                // =============================================
                // set & get methods
                // =============================================

                constexpr void set_seed(const unsigned int * s, const unsigned int& p1, const unsigned int& p2) {
                    m1 = 502;
                    m2 = 1521;
                    m3 = 4071;
                    m4 = 2107;
                    l1 = s[0] % 4096;
                    l2 = s[1] % 4096;
                    l3 = s[2] % 4096;
                    l4 = s[3] % 4096;
                    l4 = 2 * (l4 / 2) + 1;
                    n1 = 0;
                    n2 = 0;
                    n3 = p1;
                    n4 = p2;
                }

                void set_up() {
                    unsigned int seed[4];
                    unsigned int p1, p2;
                    std::ifstream file_in("primes");
                    if (file_in.is_open()) { file_in >> p1 >> p2 ; } 
                    else std::cerr << "PROBLEM: Unable to open Primes\n";
                    file_in.close();

                    file_in.open("seed.in");
                    std::string property;
                    if (file_in.is_open()) {
                        while (!file_in.eof()) {
                            file_in >> property;
                            if (property == "RANDOMSEED") {
                                file_in >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                                set_seed(seed, p1, p2);
                            }
                        }
                        file_in.close();
                    } else std::cerr << "PROBLEM: Unable to open seed.in\n";
                }

                void save_seed() const { 
                    std::ofstream file_out("seed.out");
                    if (file_out.is_open()) file_out << l1 << " " << l2 << " " << l3 << " " << l4 << "\n";
                    else std::cerr << "PROBLEM: Unable to open random.out\n";
                    file_out.close();
                }


                // =============================================
                // distributions methods
                // =============================================

                constexpr double rannyu() {
                    const double twom12{0.000244140625};
                    int i1{}, i2{}, i3{}, i4{};
                    i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
                    i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
                    i3 = l3 * m4 + l4 * m3 + n3;
                    i4 = l4 * m4 + n4;
                    l4 = i4 % 4096;
                    i3 = i3 + i4 / 4096;
                    l3 = i3 % 4096;
                    i2 = i2 + i3 / 4096;
                    l2 = i2 % 4096;
                    l1 = (i1 + i2 / 4096) % 4096;
                    return twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * l4)));
                }

                constexpr double unif(const double& min, const double& max) {
                    return min + std::fabs(max - min) * rannyu(); 
                }

                constexpr double exp(const double& lambda) {
                    return - std::log(1 - rannyu()) / lambda; 
                }

                constexpr double lorentzian(const double& mean, const double& gamma) {
                    double y = rannyu();
                    return mean + gamma * std::tan(constants::pi * (y - 0.5));
                }

                constexpr double gauss_box_muller(const double& mean, const double& sigma) {
                    double s{rannyu()}, t{rannyu()}; 
                    double x = algebra::sqrt(-2 * std::log(s)) * std::cos(2 * constants::pi * t);
                    return mean + sigma * x;
                }

                constexpr double gauss_accept_reject(const double& mean, const double& sigma) {
                    double x{}, y{}, g{}; 
                    while (true) {
                        x = unif(-5., 5.); 
                        y = rannyu(); 
                        g = exp(algebra::square(x) / 2); 
                        if (y <= g) break;
                    }
                    return mean + x * sigma;
                }


        }; // class tools::random_generator


        // class integral for evaluating the integrals
        class integral {

            private: 

                // =============================================
                // class members
                // =============================================

                double a_{}, b_{}, h_{};
                
                unsigned int steps_{};

                int sign_{}; 

                double sum_{}, integral_{}, old_integral_{}, error_{};  

                random_generator rg_;


                // =============================================
                // set methods
                // =============================================

                constexpr void steps(const unsigned int& n) { 
                    steps_ = n; 
                    h_ = std::fabs(b_ - a_) / steps_;
                }        

                constexpr void check_range() { sign_ = (a_ == b_ ? 0. : (b_ > a_ ? 1. : -1)); }
                
                constexpr void sum(const double& sum) { sum_ = sum; }

                constexpr void reset_integral() { integral_ = 0; }    

                constexpr void begin_integration(const double& a, const double& b, unsigned int n = 1000, const double& sum0 = 0) {
                    a_ = a; 
                    b_ = b; 
                    check_range(); 
                    steps(n);
                    reset_integral(); 
                    sum(sum0); 
                }


            public: 

                // =============================================
                // constructors
                // =============================================

                integral() noexcept {}
                
                ~integral() {} 

            
                // =============================================
                // get methods
                // =============================================

                constexpr double a() const { return a_; }

                constexpr double b() const { return b_; }

                constexpr int sign() const { return sign_; }

                constexpr unsigned int steps() const { return steps_; }

                constexpr double h() const { return h_; }

                constexpr double sum() const { return sum_; }

                constexpr double value() const { return integral_; }

                constexpr double error() const { return error_; }


                // =============================================
                // print methods
                // =============================================

                void print_value(const double& precision = 1.e-6) {
                    std::cout << "integral of f(x) in [" << a_ << ", " << b_ << "] = " << std::setprecision(precision) << integral_ << "\n";
                }

                void print_error(const double& precision = 1.e-6) {
                    std::cout << "error = " << std::setprecision(precision) << error_ << "\n";
                }        

                void print_integral(const double& precision = 1.e-6) {
                    print_value(precision); 
                    print_error(precision); 
                }
                
                
                // =============================================
                // integration methods
                // =============================================

                constexpr void midpoint(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    for (unsigned int i{}; i < steps_; i++) { sum_ += (f.eval(a_ + (i + 0.5) * h_)); }
                    integral_ = sum_ * h_; 
                }

                constexpr void midpoint_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    double old_integral_2{}, old_integral_3{};
                    begin_integration(a, b, 1); 
                    while (true) {
                        old_integral_3 = old_integral_2;
                        old_integral_2 = old_integral_;
                        old_integral_ = integral_; 
                        midpoint(a_, b_, f, steps_ * 2);
                        error_ = 64 * std::fabs(64 * integral_ - 84 * old_integral_ + 21 * old_integral_2 - old_integral_3) / 2835; // errore al sesto ordine
                        if (error_ < prec) break;
                    }  
                    integral_ = (4096 * integral_ - 1344 * old_integral_ + 84 * old_integral_2 - old_integral_3) / 2835; // integrale all'ottavo ordine
                }
                
                constexpr void trapexoid(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 2.);
                    for (unsigned int i{1}; i < steps_; i++) sum_ += f.eval(a_ + i * h_); 
                    integral_ = sum_ * h_; 
                } 

                constexpr void trapexoid_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    double old_integral_2{}, old_integral_3{};
                    begin_integration(a, b, 2, f.eval(a) + f.eval(b) / 2. + f.eval((a + b) / 2.)); 
                    while (true) {
                        old_integral_3 = old_integral_2;
                        old_integral_2 = old_integral_;
                        old_integral_ = integral_; 
                        trapexoid(a_, b_, f, steps_ * 2);
                        error_ = 64 * std::fabs(64 * integral_ - 84 * old_integral_ + 21 * old_integral_2 - old_integral_3) / 2835; // errore al sesto ordine
                        if (error_ < prec) break;
                    }
                    integral_ = (4096 * integral_ - 1344 * old_integral_ + 84 * old_integral_2 - old_integral_3) / 2835; // integrale all'ottavo ordine
                }

                constexpr void simpson(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    if (n % 2 == 0) begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 3.);
                    else begin_integration(a, b, n + 1);  
                    for (unsigned int i{1}; i < steps_; i++) sum_ += 2 * (1 + i % 2) * (f.eval(a_ + i * h_)) / 3.; 
                    integral_ = sum_ * h_; 
                }

                constexpr void simpson_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    double old_integral_2{}, old_integral_3{};
                    begin_integration(a, b, 2, (f.eval(a) + f.eval(b)) / 3.); 
                    while (true) {
                        old_integral_3 = old_integral_2;
                        old_integral_2 = old_integral_;
                        old_integral_ = integral_; 
                        simpson(a_, b_, f, steps_ * 2);
                        error_ = 256 * std::fabs(1024 * integral_ - 1104 * old_integral_ + 81 * old_integral_2 - old_integral_3) / 240975;
                        if (error_ < prec) break; 
                    }
                    integral_ = (1024 * integral_ - 80 * old_integral_ + old_integral_2) / 945; 
                }

                constexpr void mean(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    for (unsigned int i{}; i < n; i ++) sum_ += f.eval(rg_.unif(a, b)); 
                    integral_ = (b_ - a_) * sum_ / steps_; 
                }

                void mean_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        mean(a, b, f);
                        k.emplace_back(integral_); 
                    }
                    double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                    mean(a, b, f, static_cast<unsigned int>(algebra::square(k_mean / prec))); 
                }
        
                constexpr void hit_or_miss(const double& a, const double& b, const functions::function_base& f, const double& fmax, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    double x{}, y{}; 
                    unsigned int hits{};
                    for (unsigned int i{}; i < n; i ++) {
                        x = rg_.unif(a, b); 
                        y = rg_.unif(0., fmax);  
                        if (y <= f.eval(x)) hits++; 
                    }
                    integral_ = (b_ - a_) * fmax * hits / n; 
                }

                void hit_or_miss_fixed(const double& a, const double& b, const functions::function_base& f, const double& fmax, const double& prec = 1.e-6) {
                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        hit_or_miss(a, b, f, fmax);
                        k.emplace_back(integral_); 
                    }
                    double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                    unsigned int N = static_cast<unsigned int>(algebra::square(k_mean / prec)); 
                    hit_or_miss(a, b, f, fmax, N); 
                }


        }; // class integral


    } // namespace tools


} // namespace math


namespace physics {


    // namespace defining the units of measurements
    namespace units {
        
        // number of bits used for encoding base unit exponents 
        namespace bitwidth {

            constexpr uint32_t base_size = sizeof(uint32_t) == 8 ? 8 : 4;
            constexpr uint32_t metre{(base_size == 8) ? 8 : 4};
            constexpr uint32_t second{(base_size == 8) ? 8 : 4};
            constexpr uint32_t kilogram{(base_size == 8) ? 6 : 3};
            constexpr uint32_t ampere{(base_size == 8) ? 6 : 3};
            constexpr uint32_t candela{(base_size == 8) ? 4 : 2};
            constexpr uint32_t kelvin{(base_size == 8) ? 6 : 3};
            constexpr uint32_t mole{(base_size == 8) ? 4 : 2};

        } // namespace bitwith


        // class representing the seven SI base units 
        class unit_data {

            private: 

                signed int metre_ : bitwidth::metre;
                signed int second_ : bitwidth::second;  
                signed int kilogram_ : bitwidth::kilogram;
                signed int ampere_ : bitwidth::ampere;
                signed int kelvin_ : bitwidth::kelvin;
                signed int mole_ : bitwidth::mole;
                signed int candela_ : bitwidth::candela;  
    
                static constexpr uint32_t bits[7] = { 
                    bitwidth::metre, bitwidth::second, 
                    bitwidth::kilogram, bitwidth::ampere, 
                    bitwidth::kelvin, bitwidth::mole, 
                    bitwidth::candela
                };  


            public:

                // the seven SI base units 
                enum base {
                    Metre = 0,
                    Second = 1,
                    Kilogram = 2,
                    Ampere = 3,
                    Kelvin = 4,
                    Mole = 5,
                    Candela = 6
                };


                // =============================================
                // constructors
                // ============================================= 

                constexpr unit_data() : metre_(0), second_(0), kilogram_(0), ampere_(0), kelvin_(0), mole_(0), candela_(0) {};
                
                // constructor from powers
                constexpr unit_data(const int& metres, const int& seconds, const int& kilograms, const int& amperes, const int& kelvins,  const int& moles, const int& candelas) :
                    metre_(metres), second_(seconds), kilogram_(kilograms),
                    ampere_(amperes), kelvin_(kelvins), mole_(moles), candela_(candelas) {};

                explicit constexpr unit_data(std::nullptr_t) :
                    metre_(math::constants::max_neg(bitwidth::metre)), 
                    second_(math::constants::max_neg(bitwidth::second)), 
                    kilogram_(math::constants::max_neg(bitwidth::kilogram)), 
                    ampere_(math::constants::max_neg(bitwidth::ampere)),
                    kelvin_(math::constants::max_neg(bitwidth::kelvin)), 
                    mole_(math::constants::max_neg(bitwidth::mole)), 
                    candela_(math::constants::max_neg(bitwidth::candela)) {}


                // =============================================
                // operators
                // ============================================= 

                // perform a multiply operation by adding the powers together
                constexpr unit_data operator*(const unit_data& other) const {
                    return { 
                        metre_ + other.metre_,
                        second_ + other.second_,
                        kilogram_ + other.kilogram_,
                        ampere_ + other.ampere_,
                        kelvin_ + other.kelvin_,
                        mole_ + other.mole_,
                        candela_ + other.candela_,
                    };
                }

                // perform a division operation by subtract the powers together
                constexpr unit_data operator/(const unit_data& other) const {
                    return { 
                        metre_ - other.metre_,
                        second_ - other.second_,
                        kilogram_ - other.kilogram_,
                        ampere_ - other.ampere_,
                        kelvin_ - other.kelvin_,
                        mole_ - other.mole_,
                        candela_ - other.candela_,
                    };
                }

                // comparison operators
                constexpr bool operator==(const unit_data& other) const {
                    return metre_ == other.metre_ && second_ == other.second_ &&
                        kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                        candela_ == other.candela_ && kelvin_ == other.kelvin_ && mole_ == other.mole_;         
                }

                // definitely not a comparison operators
                constexpr bool operator!=(const unit_data& other) const { return !(*this == other); }


                // =============================================
                // operations
                // ============================================= 

                // invert the unit
                constexpr unit_data inv() const {
                    return { 
                        -metre_,
                        -second_,
                        -kilogram_,
                        -ampere_,
                        -kelvin_,
                        -mole_,
                        -candela_
                    };
                }

                // take a unit_data to some power
                constexpr unit_data pow(const int& power) const { 
                    return { 
                        metre_ * power,
                        second_ * power,
                        kilogram_ * power,
                        ampere_ * power,
                        kelvin_ * power,
                        mole_ * power,
                        candela_ * power
                    };
                }

                constexpr unit_data square() const { return pow(2); }

                constexpr unit_data cube() const { return pow(3); }
                
                // take some root of a unit_data
                constexpr unit_data root(const int& power) const {
                    if (has_valid_root(power)) return unit_data(metre_ / power,
                                                                (second_ / power) + root_hertz_modifier(power),
                                                                kilogram_ / power,
                                                                ampere_ / power,
                                                                kelvin_ / power,
                                                                mole_ / power,
                                                                candela_ / power);
                    else exit(-11);
                }

                constexpr unit_data sqrt() const { return root(2); }

                constexpr unit_data cbrt() const { return root(3); }
                                    

                // =============================================
                // get methods
                // =============================================
                
                // check if the units have the same base unit 
                constexpr bool has_same_base(const unit_data& other) const { return *this == other; }

                constexpr int metre() const { return metre_; }

                constexpr int second() const { return second_; }

                constexpr int kg() const { return kilogram_; }
                                    
                constexpr int ampere() const { return ampere_; }
                
                constexpr int kelvin() const { return kelvin_; }
                
                constexpr int mole() const { return mole_; }
                
                constexpr int candela() const { return candela_; }

                constexpr unit_data data() const { return *this; }

                constexpr void print_data() const {
                    if (metre_ != 0 && metre_ != 1) std::cout << "m^" << metre_; 
                    else if (metre_ == 1) std::cout << "m";
                    if (second_ != 0 && second_ != 1) std::cout << "s^" << second_; 
                    else if (second_ == 1) std::cout << "s"; 
                    if (kilogram_ != 0 && kilogram_ != 1) std::cout << "kg^" << kilogram_; 
                    else if (kilogram_ == 1) std::cout << "kg"; 
                    if (ampere_ != 0 && ampere_ != 1) std::cout << "A^" << ampere_; 
                    else if (ampere_ == 1) std::cout << "A"; 
                    if (kelvin_ != 0 && kelvin_ != 1) std::cout << "K^" << kelvin_; 
                    else if (kelvin_ == 1) std::cout << "K";
                    if (mole_ != 0 && mole_ != 1) std::cout << "mol^" << mole_; 
                    else if (mole_ == 1) std::cout << "mol"; 
                    if (candela_ != 0 && candela_ != 1) std::cout << "cd^" << candela_; 
                    else if (candela_ == 1) std::cout << "cd"; 
                }


            private: 

                // check if the base_unit has a valid root
                constexpr bool has_valid_root(const int& power) const {
                    return metre_ % power == 0 && second_ % power == 0 &&
                        kilogram_ % power == 0 && ampere_ % power == 0 &&
                        candela_ % power == 0 && kelvin_ % power == 0 && mole_ % power == 0;
                }      
                
                // to handle a few weird operations that operate on square_root Hz
                constexpr int root_hertz_modifier(const int& power) const {
                    return (second_ * power == 0 || power % 2 != 0) ? 0 : (power / 2) * ((second_ < 0) || (power < 0) ? 9 : -9);
                }    


        }; // class unit_base


        // class defining a basic unit module 
        class unit_prefix {

            protected:  

                // =============================================
                // class members
                // ============================================= 
                
                double multiplier_{};  


            public:  

                // =============================================
                // constructors
                // ============================================= 
                                    
                constexpr unit_prefix(const double& mult = 1.0) noexcept : multiplier_{mult} {} 


                // =============================================
                // operators
                // ============================================= 

                // perform a multiply operation by adding the powers together
                constexpr unit_prefix operator*(const unit_prefix& other) const { return { multiplier_ * other.multiplier_ }; }

                // perform a division operation by subtract the powers together
                constexpr unit_prefix operator/(const unit_prefix& other) const { return { multiplier_ / other.multiplier_ }; }

                // comparison operators
                constexpr bool operator==(const unit_prefix& other) const {
                    if (multiplier_ == other.multiplier_) return true;    
                    else return math::tools::compare_round_equals(multiplier_, other.multiplier_);
                }

                // definitely not a comparison operators
                constexpr bool operator!=(const unit_prefix& other) const { return !(*this == other); }


                // =============================================
                // operations
                // ============================================= 

                // invert the unit_prefix
                constexpr unit_prefix inv() const { return { 1 / multiplier_ }; }

                // take a unit_prefix to some power
                constexpr unit_prefix pow(const int& power) const { return { std::pow(multiplier_, power) }; }
                
                constexpr unit_prefix square() const { return { std::pow(multiplier_, 2) }; }
                
                constexpr unit_prefix cube() const { return { std::pow(multiplier_, 3) }; }

                // take some root of a unit_prefix
                constexpr unit_prefix root(const int& power) const { return { math::algebra::root(multiplier_, power) }; }

                constexpr unit_prefix sqrt() const { return { math::algebra::root(multiplier_, 2) }; }
                
                constexpr unit_prefix cbrt() const { return { math::algebra::root(multiplier_, 3) }; }
                
                
                // =============================================
                // check methods
                // ============================================= 

                // check if the units have the same prefix unit 
                constexpr bool has_same_prefix(const unit_prefix& other) const { return *this == other; }

                // check if the multiplier is nan
                constexpr bool is_nan(const unit_prefix& pref) { return std::isnan(pref.multiplier_); }

                // checks that the multiplier is finite
                constexpr bool is_finite(const unit_prefix& pref) { return std::isfinite(pref.multiplier_); }

                // check if the multiplier is infinite
                constexpr bool is_inf(const unit_prefix& pref) { return std::isinf(pref.multiplier_); }


                // =============================================
                // get methods
                // ============================================= 

                // get the multiplier
                constexpr double multiplier() const { return multiplier_; }

                // get the prefix
                constexpr unit_prefix prefix() const { return *this; }

                // print the prefix
                constexpr void print_prefix() const { 
                    if (multiplier_ == 1e-1) std::cout << 'd'; 
                    if (multiplier_ == 1e-2) std::cout << 'c'; 
                    if (multiplier_ == 1e-3) std::cout << 'm'; 
                    if (multiplier_ == 1e-6) std::cout << 'u'; 
                    if (multiplier_ == 1e-9) std::cout << 'n'; 
                    if (multiplier_ == 1e-12) std::cout << 'p'; 
                    if (multiplier_ == 1e-15) std::cout << 'f'; 
                    if (multiplier_ == 1e-18) std::cout << 'a'; 
                    if (multiplier_ == 1e-21) std::cout << 'z'; 
                    if (multiplier_ == 1e-24) std::cout << 'y'; 
                    if (multiplier_ == 1e2) std::cout << 'h'; 
                    if (multiplier_ == 1e3) std::cout << 'k'; 
                    if (multiplier_ == 1e6) std::cout << 'M'; 
                    if (multiplier_ == 1e9) std::cout << 'G'; 
                    if (multiplier_ == 1e12) std::cout << 'T'; 
                    if (multiplier_ == 1e15) std::cout << 'P'; 
                    if (multiplier_ == 1e18) std::cout << 'E'; 
                    if (multiplier_ == 1e21) std::cout << 'Z'; 
                    if (multiplier_ == 1e24) std::cout << 'Y';
                }


        }; // class unit_prefix


        // class defining a basic unit module 
        class unit : public unit_data, public unit_prefix {

            public:

                // =============================================
                // constructors & destructor
                // ============================================= 

                // constructor from unit_data 
                explicit constexpr unit(const unit_data& data = unit_data(), const unit_prefix& prefix = unit_prefix()) noexcept : unit_data(data), unit_prefix(prefix) {}

                // constructor from unit
                constexpr unit(const unit& other) noexcept : unit_data(other.data()), unit_prefix(other.prefix()) {} 

                // constructor from unit and unit_prefix
                constexpr unit(const unit& other, const unit_prefix& prefix) noexcept : unit_data(other.data()), unit_prefix(prefix) {} 

                // constructor from unit_prefix and unit
                constexpr unit(const unit_prefix& prefix, const unit& other) noexcept : unit_data(other.data()), unit_prefix(prefix) {} 

                // default destructor
                ~unit() = default; 


                // =============================================
                // operators
                // ============================================= 

                constexpr unit operator*(const unit& other) const { return unit(unit_data::operator*(other.data()), unit_prefix::operator*(other.prefix())); }

                constexpr unit operator/(const unit& other) const { return unit(unit_data::operator/(other.data()), unit_prefix::operator/(other.prefix())); }                    

                constexpr bool operator==(const unit& other) const {
                    if (unit_data::operator!=(other.data())) { return false; }
                    if (unit_prefix::operator==(other.prefix())) { return true; }
                    return math::tools::compare_round_equals(multiplier(), other.multiplier());
                }

                constexpr bool operator!=(const unit& other) const { return !operator==(other); }


                // =============================================
                // operations
                // ============================================= 

                constexpr unit inv() const { return unit(unit_data::inv(), unit_prefix::inv()); }

                constexpr unit pow(const int& power) const { return unit(unit_data::pow(power), unit_prefix::pow(power)); }

                constexpr unit square() const { return unit(unit_data::square(), unit_prefix::square()); }

                constexpr unit cube() const { return unit(unit_data::cube(), unit_prefix::cube()); }

                constexpr unit root(const int& power) const { return unit(unit_data::root(power), unit_prefix::root(power)); }

                constexpr unit sqrt() const { return unit(unit_data::sqrt(), unit_prefix::sqrt()); }

                constexpr unit cbrt() const { return unit(unit_data::cbrt(), unit_prefix::cbrt()); }


                // =============================================
                // get methods
                // ============================================= 

                // test for exact numerical equivalence
                constexpr bool is_exactly_the_same(const unit& other) const { 
                    return unit_data::operator==(other.data()) && unit_prefix::operator==(other.prefix());
                }

                // get the unit
                constexpr unit as_unit() const { return *this; }

                // print the unit
                constexpr void print_unit() const {
                    print_prefix(); 
                    print_data();
                }


        }; // class unit     


        // base
        namespace base {

            constexpr unit_data metre(1, 0, 0, 0, 0, 0, 0);
            constexpr unit_data second(0, 1, 0, 0, 0, 0, 0);
            constexpr unit_data kilogram(0, 0, 1, 0, 0, 0, 0);
            constexpr unit_data Ampere(0, 0, 0, 1, 0, 0, 0);
            constexpr unit_data Kelvin(0, 0, 0, 0, 1, 0, 0);
            constexpr unit_data mol(0, 0, 0, 0, 0, 1, 0);
            constexpr unit_data candela(0, 0, 0, 0, 0, 0, 1);

        } // namespace SI


        // prefix
        namespace prefix {

            constexpr unit_prefix deci(1e-1); 
            constexpr unit_prefix centi(1e-2); 
            constexpr unit_prefix milli(1e-3); 
            constexpr unit_prefix micro(1e-6); 
            constexpr unit_prefix nano(1e-9); 
            constexpr unit_prefix pico(1e-12); 
            constexpr unit_prefix femto(1e-15); 
            constexpr unit_prefix atto(1e-18); 
            constexpr unit_prefix zepto(1e-21); 
            constexpr unit_prefix yocto(1e-24); 
            constexpr unit_prefix hecto(1e2); 
            constexpr unit_prefix kilo(1e3); 
            constexpr unit_prefix mega(1e6); 
            constexpr unit_prefix giga(1e9); 
            constexpr unit_prefix tera(1e12); 
            constexpr unit_prefix peta(1e15); 
            constexpr unit_prefix exa(1e18); 
            constexpr unit_prefix zetta(1e21); 
            constexpr unit_prefix yotta(1e24); 

        } // namespace prefix


        // special units 
        namespace special {

            // some unitless numbers
            unit one;
            unit hundred = unit(100.0, one);
            unit ten = unit(10.0, one);
            unit percent(one, 0.01);
            unit infinite(unit_data(0, 0, 0, 0, 0, 0, 0), math::constants::infinity);
            unit neginfinite(unit_data(0, 0, 0, 0, 0, 0, 0), -math::constants::infinity);
            unit nan(unit_data(0, 0, 0, 0, 0, 0, 0), math::constants::invalid_conversion);

            // some specialized units
            unit defunit(unit_data(0, 0, 0, 0, 0, 0, 0));
            unit invalid(unit_data(nullptr), math::constants::invalid_conversion);
            unit error(unit_data(nullptr));                
        
        } // namespace special


        // SI units
        namespace SI {
                        
            // SI units
            unit m = unit(base::metre);
            unit s = unit(base::second);
            unit kg = unit(base::kilogram);
            unit A = unit(base::Ampere);
            unit K = unit(base::Kelvin);
            unit mol = unit(base::mol); 
            unit cd = unit(base::candela);

            // distance units
            unit km(prefix::kilo, m);
            unit dm(prefix::deci, m);
            unit cm(prefix::centi, m);
            unit mm(prefix::milli, m);
            unit um(prefix::micro, m);
            unit nm(prefix::nano, m);

            // time units
            unit ms(prefix::milli, s);
            unit us(prefix::micro, s);
            unit ns(prefix::nano, s);
            unit min(60.0, s);
            unit hr(60.0, min);
            unit day(24.0, hr);
            
            unit radians = unit(180 / math::constants::pi, special::one); 
            unit rad = radians;

            unit mps(m / s);
            unit mpss(m / s.pow(2)); 
            unit kmps(km / s);
            unit kmpss(km / s.pow(2)); 


            // derived_units
            namespace derived {

                unit hertz(SI::s.inv());
                unit Hz = hertz;

                unit volt(unit_data(2, -3, 1, -1, 0, 0, 0));
                unit V = volt;

                unit newton(SI::kg * SI::mpss);
                unit N = newton;

                unit Pa(unit_data(-1, -2, 1, 0, 0, 0, 0));
                unit pascal = Pa;

                unit joule(unit_data(2, -2, 1, 0, 0, 0, 0));
                unit J = joule;

                unit watt(unit_data(2, -3, 1, 0, 0, 0, 0));
                unit W = watt;

                unit coulomb(unit_data(0, 1, 0, 1, 0, 0, 0));
                unit C = coulomb;

                unit farad(unit_data(-2, 4, -1, 2, 0, 0, 0));
                unit F = farad;

                unit weber(unit_data(2, -2, 1, -1, 0, 0, 0));
                unit Wb = weber;

                unit tesla(unit_data(0, -2, 1, -1, 0, 0, 0));
                unit T = tesla;

                unit henry(unit_data(2, -2, 1, -2, 0, 0, 0));                    
                unit H = henry;


            } // namespace derived


        } // namespace SI


        // unit yr(8760.0, hr, "yr");  // median calendar year;
        // unit sday(365.24 / 366.24, day, "sday");  // sidereal day
        // unit syr(365.256363004, day, "syr");  // sidereal year

        // // mass units
        // unit g(prefix::milli, kg, "g");
        // unit mg(prefix::micro, kg, "mg");

        // // volume units
        // unit L{1, dm.pow(3), "L"};
        // unit dL{0.1, L, "dL"};
        // unit cL{0.01, L, "cL"};
        // unit mL{0.001, L, "mL"};                    


    } // namespace units


} // namespace physics


namespace math {


    // namespace defining some usefull operations
    namespace tools {

        // generate a conversion factor between two units in a constexpr function, the
        // units will only convert if they have the same base unit
        template<typename UX, typename UX2>
        constexpr double quick_convert(UX start, UX2 result) { return quick_convert(1.0, start, result); }

        // generate a conversion factor between two units in a constexpr function, the units will only convert if they have the same base unit
        template<typename UX, typename UX2>
        constexpr double quick_convert(double val, const UX& start, const UX2& result) {
            static_assert(std::is_same<UX, physics::units::unit>::value || std::is_same<UX, physics::units::unit>::value, "convert argument types must be unit");
            static_assert(std::is_same<UX2, physics::units::unit>::value || std::is_same<UX2, physics::units::unit>::value, "convert argument types must be unit");
            return (start.data() == result.data()) ?  val * start.multiplier() / result.multiplier() : constants::invalid_conversion;
        }

        // convert a value from one unit base to another
        template<typename UX, typename UX2>
        constexpr double convert(double val, const UX& start, const UX2& result) {
            static_assert(std::is_same<UX, physics::units::unit>::value || std::is_same<UX, physics::units::unit>::value, "convert argument types must be unit");
            static_assert(std::is_same<UX2, physics::units::unit>::value || std::is_same<UX2, physics::units::unit>::value, "convert argument types must be unit");     
            if (start == result) { return val; }
            if (start.data() == result.data()) { return val * start.multiplier() / result.multiplier(); }
            auto base_start = start.data();
            auto base_result = result.data();
            if (base_start.has_same_base(base_result)) { return val * start.multiplier() / result.multiplier(); }
            if (base_start.has_same_base(base_result.inv())) { return 1.0 / (val * start.multiplier() * result.multiplier()); }
            return constants::invalid_conversion;
        }     
                        
        // generate a conversion factor between two physics::units::units
        template<typename UX, typename UX2>
        constexpr double convert(const UX& start, const UX2& result) { return convert(1.0, start, result); }


    } // namespace tools


} // namespace math


namespace physics {


    // namespace defining some usefull measurement structures
    namespace measurements {
        
        
        using namespace units::SI;
        using namespace units::SI::derived;


        // measurement class with a value and an unit
        class measurement {

            protected:

                // =============================================                                                                                         
                // class members
                // =============================================  

                double value_{};

                units::unit unit_;


            public:

                // =============================================                                                                                         
                // constructors & destructor 
                // =============================================  

                // default constructor
                constexpr measurement() noexcept {};

                // constructor from a value and a unit 
                explicit constexpr measurement(const double& value, const units::unit& unit) noexcept : value_{value}, unit_{unit} {}

                // constructor from a measurement
                constexpr measurement(const measurement& other) : value_{other.value_}, unit_{other.unit_} {}
                    
                // default destructor
                ~measurement() = default;


                // =============================================                                                                                         
                // operators
                // =============================================  

                constexpr measurement operator+(const measurement& other) const { return measurement(value_ + other.value_as(unit_), unit_); }
                
                constexpr measurement operator+(const double& val) const { return measurement(value_ + val, unit_); }

                constexpr measurement operator-() const { return measurement(- value_, unit_); }

                constexpr measurement operator-(const measurement& other) const { return measurement(value_ - other.value_as(unit_), unit_); }

                constexpr measurement operator-(const double& val) const { return measurement(value_ - val, unit_); }

                constexpr measurement operator*(const double& val) const { return measurement(value_ * val, unit_); }

                constexpr measurement operator*(const measurement& other) const { return measurement(value_ * other.value_, unit_ * other.unit_); }
                                
                constexpr measurement operator/(const double& val) const { return measurement(value_ / val, unit_); } 

                constexpr measurement operator/(const measurement& other) const { return measurement(value_ / other.value_, unit_ / other.unit_); }

                constexpr measurement& operator=(const double& val) noexcept {
                    value_ = val;
                    return *this;
                }

                constexpr measurement& operator=(const measurement& other) noexcept {
                    value_ = other.value_;
                    unit_ = other.unit_;
                    return *this;
                }     

                constexpr bool operator==(const double& val) const { return (value_ == val) ? true : math::tools::compare_round_equals(value_, val); }

                constexpr bool operator==(const measurement& other) const { return value_equality_check((unit_ == other.units()) ? other.value() : other.value_as(unit_)); }

                constexpr bool operator!=(const double& val) const { return !operator==(val); }

                constexpr bool operator!=(const measurement& other) const { return !value_equality_check((unit_ == other.units()) ? other.value() : other.value_as(unit_)); }

                constexpr bool operator>(const double& val) const { return value_ > val; }

                constexpr bool operator>(const measurement& other) const { return value_ > other.value_as(unit_); }

                constexpr bool operator<(const double& val) const { return value_ < val; }

                constexpr bool operator<(const measurement& other) const { return value_ < other.value_as(unit_); }

                constexpr bool operator>=(const double& val) const { return (value_ >= val) ? true : operator==(val); }

                constexpr bool operator>=(const measurement& other) const { return (value_ > other.value_as(unit_)) ? true : value_equality_check(other.value_as(unit_)); }

                constexpr bool operator<=(const double& val) const { return value_ <= val ? true : operator==(val); }

                constexpr bool operator<=(const measurement& other) const { return (value_ < other.value_as(unit_)) ? true : value_equality_check(other.value_as(unit_)); }       


                // =============================================
                // operations
                // ============================================= 

                constexpr measurement inv() const { return measurement(1 / value_, unit_.inv()); }

                constexpr measurement pow(const int& power) const { return measurement(std::pow(value_, power), unit_.pow(power)); }

                constexpr measurement square() const { return measurement(math::algebra::square(value_), unit_.square()); }

                constexpr measurement cube() const { return measurement(math::algebra::cube(value_), unit_.cube()); }

                constexpr measurement root(const int& power) const { return measurement(math::algebra::root(value_, power), unit_.root(power)); }

                constexpr measurement sqrt() const { return measurement(math::algebra::sqrt(value_), unit_.sqrt()); }

                constexpr measurement cbrt() const { return measurement(math::algebra::cbrt(value_), unit_.cbrt()); }

                
                // =============================================                                                                                         
                // set & get, convert & print methods
                // =============================================  

                // set the numerical value of the measurement
                constexpr void value(const double& value) { value_ = value; }

                // get the numerical value of the measurement
                constexpr double value() const { return value_; }                        
                
                // get the numerical value as a particular unit type
                constexpr double value_as(const units::unit& desired_units) const { return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units); }

                // get the unit component of a measurement
                constexpr units::unit units() const { return unit_; }

                // convert a unit to have a new base
                constexpr measurement convert_to(const units::unit& newUnits) const { return measurement(math::tools::convert(value_, unit_, newUnits), newUnits); }

                // print the measurement
                void print_measurement() const { 
                    // if (unit_.multiplier()!= 1) std::cout << std::setprecision(unit_.multiplier()) << value_ << " "; 
                    std::cout << value_ << " ";
                    unit_.print_unit(); 
                    std::cout << "\n"; 
                }

            private:

                // does a numerical equality check on the value accounting for tolerances
                constexpr bool value_equality_check(const double& val) const { return (value_ == val) ? true : math::tools::compare_round_equals(value_, val); } 


        }; // class measurement


        // measurement class with an uncertain value and an unit
        class uncertain_measurement {

            private:

                // =============================================                                                                                         
                // class members
                // =============================================  

                double value_{}, uncertainty_{};

                units::unit unit_;       
            

            public: 

                // =============================================                                                                                         
                // constructors & destructor 
                // =============================================  

                // default constructor
                constexpr uncertain_measurement() noexcept {};

                // constructor from a value, uncertainty, and unit
                explicit constexpr uncertain_measurement(const double& val, const double& uncertainty_val, const units::unit& unit) noexcept :
                    value_(val), uncertainty_(uncertainty_val), unit_(unit) {}

                // constructpr from a value and an unit, assuming the uncertainty is 0
                explicit constexpr uncertain_measurement(const double& val, const units::unit& unit) noexcept :
                    value_(val), unit_(unit) {}

                // constructor from a measurement and uncertainty value
                explicit constexpr uncertain_measurement(const measurement& other, const double& uncertainty_val) noexcept : 
                    value_(other.value()), uncertainty_(uncertainty_val), unit_(other.units()) {}

                // constructor from a uncertain_measurement
                constexpr uncertain_measurement(const uncertain_measurement& other) :
                    value_(other.value_), uncertainty_(other.uncertainty_), unit_(other.unit_) {}

                // default constructor
                ~uncertain_measurement() = default;
                

                // =============================================                                                                                         
                // operators
                // =============================================  

                // compute a product and calculate the new uncertainties using the root sum of squares(rss) method
                constexpr uncertain_measurement operator*(const uncertain_measurement& other) const {
                    double tval1 = uncertainty_ / value_;
                    double tval2 = other.uncertainty_ / other.value_;
                    double ntol = math::algebra::sqrt(tval1 * tval1 + tval2 * tval2);
                    double nval = value_ * other.value_;
                    return uncertain_measurement(nval, nval * ntol, unit_ * other.units());
                }

                // perform a multiplication with uncertain measurements using the simple method for uncertainty propagation
                constexpr uncertain_measurement simple_product(const uncertain_measurement& other) const {
                    double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                    double nval = value_ * other.value_;
                    return uncertain_measurement(nval, nval * ntol, unit_ * other.units());
                }

                // multiply with another measurement equivalent to uncertain_measurement multiplication with 0 uncertainty
                constexpr uncertain_measurement operator*(const measurement& other) const {
                    return uncertain_measurement(value() * other.value(), other.value() * uncertainty(), unit_ * other.units());
                }

                constexpr uncertain_measurement operator*(const double& val) const { return uncertain_measurement(value_ * val, uncertainty_ * val, unit_); }

                // compute a unit division and calculate the new uncertainties using the root sum of squares(rss) method
                constexpr uncertain_measurement operator/(const uncertain_measurement& other) const {
                    double tval1 = uncertainty_ / value_;
                    double tval2 = other.uncertainty_ / other.value_;
                    double ntol = math::algebra::sqrt(tval1 * tval1 + tval2 * tval2);
                    double nval = value_ / other.value_;
                    return uncertain_measurement(nval, nval * ntol, unit_ / other.units());
                }

                // division operator propagate uncertainty using simple method
                constexpr uncertain_measurement simple_divide(const uncertain_measurement& other) const {
                    double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                    double nval = value_ / other.value_;
                    return uncertain_measurement(nval, nval * ntol, unit_ / other.units());
                }

                constexpr uncertain_measurement operator/(const measurement& other) const {
                    return uncertain_measurement(static_cast<double>(value() / other.value()), static_cast<double>(uncertainty() / other.value()), unit_ / other.units());
                }

                constexpr uncertain_measurement operator/(const double& val) const {
                    return uncertain_measurement(value_ / val, uncertainty_ / val, unit_);
                }

                // compute a unit addition and calculate the new uncertainties using the root sum of squares(rss) method
                constexpr uncertain_measurement operator+(const uncertain_measurement& other) const {
                    auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                    double ntol = math::algebra::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                    return uncertain_measurement(value_ + cval * other.value_, ntol, unit_);
                }

                // compute a unit addition and calculate the new uncertainties using the simple uncertainty summation method
                constexpr uncertain_measurement simple_add(const uncertain_measurement& other) const {
                    auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                    double ntol = uncertainty_ + other.uncertainty_ * cval;
                    return uncertain_measurement(value_ + cval * other.value_, ntol, unit_);
                }

                constexpr uncertain_measurement operator+(const measurement& other) const {
                    return uncertain_measurement(value_ + other.value_as(unit_), uncertainty_, unit_);
                }

                // compute a unit subtraction and calculate the new uncertainties using the root sum of squares(rss) method
                constexpr uncertain_measurement operator-(const uncertain_measurement& other) const {
                    auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                    double ntol = math::algebra::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                    return uncertain_measurement(value_ - cval * other.value_, ntol, unit_);
                }

                // compute a unit subtraction and calculate the new uncertainties using the simple uncertainty summation method
                constexpr uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                    auto cval = math::tools::convert(other.unit_, unit_);
                    double ntol = uncertainty_ + other.uncertainty_ * cval;
                    return uncertain_measurement(value_ - cval * other.value_, ntol, unit_);
                }

                constexpr uncertain_measurement operator-(const measurement& other) const {
                    return uncertain_measurement(value_ - other.value_as(unit_), uncertainty_, unit_);
                }

                constexpr uncertain_measurement& operator=(const uncertain_measurement& other) noexcept {
                    value_ = other.value_;
                    uncertainty_ = other.uncertainty_; 
                    unit_ = other.unit_;
                    return *this;
                }     

                constexpr bool operator==(const uncertain_measurement& other) const { return (simple_subtract(other) == measurement(0.0, unit_)); }

                constexpr bool operator==(const measurement& other) const {
                    if (uncertainty_ == 0.0F) { return (value_ == other.value_as(unit_)) ? true : math::tools::compare_round_equals(value_, other.value_as(unit_)); }
                    return (other.value_as(unit_) >= (value_ - uncertainty_) && other.value_as(unit_) <= (value_ + uncertainty_));
                }

                constexpr bool operator!=(const uncertain_measurement& other) const { return !operator==(other); }

                constexpr bool operator!=(const measurement& other) const { return !operator==(other); }

                constexpr bool operator>(const uncertain_measurement& other) const { return value_ > other.value_as(unit_); }

                constexpr bool operator>(const measurement& other) const { return value_ > other.value_as(unit_); }

                constexpr bool operator>(const double& val) const { return value_ > val; }

                constexpr bool operator<(const uncertain_measurement& other) const { return value_ < other.value_as(unit_); }

                constexpr bool operator<(const measurement& other) const { return value_ < other.value_as(unit_); }

                constexpr bool operator<(const double& val) const { return value_ < val; }

                constexpr bool operator>=(const uncertain_measurement& other) const {
                    return (simple_subtract(other).value_ >= 0.0F) ? true : (simple_subtract(other) == measurement(0.0, unit_));
                }

                constexpr bool operator>=(const measurement& other) const {
                    return (value_ >= other.value_as(unit_)) ? true : operator==(measurement(other.value_as(unit_), unit_));
                }

                constexpr bool operator>=(const double& val) const { return value_ >= val; }

                constexpr bool operator<=(const uncertain_measurement& other) const {
                    return (simple_subtract(other).value_ <= 0.0F) ? true : (simple_subtract(other) == measurement(0.0, unit_));
                }

                constexpr bool operator<=(const measurement& other) const {
                    return (value_ <= other.value_as(unit_)) ? true : operator==(measurement(other.value_as(unit_), unit_));
                }

                constexpr bool operator<=(const double& val) const { return value_ <= val; }


                // =============================================
                // operations
                // ============================================= 

                constexpr uncertain_measurement inv() const { return uncertain_measurement(1 / value_, uncertainty_, unit_.inv()); }

                constexpr uncertain_measurement pow(const int& power) const { return uncertain_measurement(std::pow(value_, power), uncertainty_ * power, unit_.pow(power)); }

                constexpr uncertain_measurement square() const { return uncertain_measurement(math::algebra::square(value_), uncertainty_ * 2, unit_.square()); }

                constexpr uncertain_measurement cube() const { return uncertain_measurement(math::algebra::cube(value_), uncertainty_ * 3, unit_.cube()); }

                constexpr uncertain_measurement root(const int& power) const { return uncertain_measurement(math::algebra::root(value_, power), uncertainty_ / power, unit_.root(power)); }

                constexpr uncertain_measurement sqrt() const { return uncertain_measurement(math::algebra::sqrt(value_), uncertainty_ / 2, unit_.sqrt()); }

                constexpr uncertain_measurement cbrt() const { return uncertain_measurement(math::algebra::cbrt(value_), uncertainty_ / 3, unit_.cbrt()); }


                // =============================================                                                                                         
                // get & set & print methods
                // =============================================  
                
                // set the uncertainty
                constexpr void uncertainty(const double& uncert) { uncertainty_ = uncert; }

                // set the uncertainty
                constexpr uncertain_measurement& uncertainty(double newUncertainty) {
                    uncertainty_ = newUncertainty;
                    return *this;
                }

                // set the uncertainty
                constexpr uncertain_measurement& uncertainty(const measurement& newUncertainty) {
                    uncertainty_ = newUncertainty.value_as(unit_);
                    return *this;
                }

                // get the uncertainty as a separate measurement
                constexpr measurement uncertainty_measurement() const { return uncertain_measurement(uncertainty(), unit_); }

                // get the numerical value 
                constexpr double value() const { return value_; }

                // get the numerical value of the uncertainty
                constexpr double uncertainty() const { return uncertainty_; }

                // get the underlying units value
                constexpr units::unit units() const { return unit_; }

                // get the numerical value as a particular unit type
                constexpr double value_as(const units::unit& desired_units) const {
                    return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units);
                }

                // get the numerical value of the uncertainty as a particular unit
                constexpr double uncertainty_as(const units::unit& desired_units) const {
                    return (unit_ == desired_units) ? uncertainty_ : math::tools::convert(uncertainty_, unit_, desired_units);
                }

                // cast operator to a measurement
                constexpr operator measurement() const { return measurement(value(), unit_); }

                // convert a unit to have a new base
                constexpr uncertain_measurement convert_to(const units::unit& newUnits) const {
                    auto cval = math::tools::convert(unit_, newUnits);
                    return uncertain_measurement(cval * value_, uncertainty_ * cval, newUnits);
                }

                // print the uncertain measurement
                void print_uncertain_measurement() const { 
                    // if (unit_.multiplier() != 1) std::cout << std::setprecision(unit_.multiplier()) << value() << " ± " << uncertainty() << " "; 
                    std::cout << value() << " ± " << uncertainty() << " "; 
                    unit_.print_unit(); 
                    std::cout << "\n";
                }


        }; // class uncertain_measurement


        constexpr measurement operator*(const double& val, const measurement& meas) { return meas * val; }
        
        constexpr measurement operator*(const int& val, const measurement& meas) { return meas * val; }

        constexpr measurement operator*(const double& val, const units::unit& unit_base) { return measurement(val, unit_base); }
        
        constexpr measurement operator*(const units::unit& unit_base, const double& val) { return measurement(val, unit_base); }

        constexpr measurement operator/(const double& val, const measurement& meas) { return measurement(val / meas.value(), meas.units().inv()); }
        
        constexpr measurement operator/(const double& val, const units::unit& unit_base) { return measurement(val, unit_base.inv()); }
        
        constexpr measurement operator/(const units::unit& unit_base, const double& val) { return measurement(1.0 / val, unit_base); }

        constexpr bool operator==(const measurement& other, const uncertain_measurement& v2) { return v2 == other; }

        constexpr bool operator!=(const measurement& other, const uncertain_measurement& v2) { return v2 != other; }

        constexpr bool operator>(const measurement& other, const uncertain_measurement& v2) { return other.value() > v2.value(); }

        constexpr bool operator<(const measurement& other, const uncertain_measurement& v2) { return other.value() < v2.value(); }

        constexpr bool operator>=(const measurement& other, const uncertain_measurement& v2) { return (other > v2) ? true : (v2 == other); }

        constexpr bool operator<=(const measurement& other, const uncertain_measurement& v2) { return (other < v2) ? true : (v2 == other); }

        constexpr uncertain_measurement operator*(const measurement& v1, const uncertain_measurement& v2) { return v2.operator*(v1); }

        constexpr uncertain_measurement operator*(const double& v1, const uncertain_measurement& v2) { return v2.operator*(v1); }

        constexpr uncertain_measurement operator/(const measurement& v1, const uncertain_measurement& v2) {
            double ntol = v2.uncertainty() / v2.value();
            double nval = v1.value() / v2.value();
            return uncertain_measurement(nval, nval * ntol, v1.units() / v2.units());
        }

        constexpr uncertain_measurement operator/(const double& v1, const uncertain_measurement& v2) {
            double ntol = v2.uncertainty() / v2.value();
            double nval = v1 / v2.value();
            return uncertain_measurement(nval, nval * ntol, v2.units().inv());
        }

        constexpr uncertain_measurement operator+(const measurement& v1, const uncertain_measurement& v2) {
            double cval = math::tools::convert(v2.units(), v1.units());
            double ntol = v2.uncertainty() * cval;
            return uncertain_measurement(v1.value() + cval * v2.value(), ntol, v1.units());
        }

        constexpr uncertain_measurement operator-(const measurement& v1, const uncertain_measurement& v2) {
            double cval = math::tools::convert(v2.units(), v1.units());
            double ntol = v2.uncertainty() * cval;
            return uncertain_measurement(v1.value() - cval * v2.value(), ntol, v1.units());
        }


    } // namespace measurements


} // namespace physics


namespace std {
    
    // operations amongs std::vector<T>
    std::vector<physics::measurements::measurement> operator+(const std::vector<physics::measurements::measurement>& v1, const std::vector<physics::measurements::measurement>& v2) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + v2[i]); 
        return vec;
    }

    std::vector<physics::measurements::measurement> operator+(const std::vector<physics::measurements::measurement>& v1, const physics::measurements::measurement& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + val);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator+(const physics::measurements::measurement& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + val);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator+(const std::vector<physics::measurements::measurement>& v1, const double& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + val);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator+(const double& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + val);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator-(const std::vector<physics::measurements::measurement>& v1, const std::vector<physics::measurements::measurement>& v2) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] - v2[i]);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator-(const std::vector<physics::measurements::measurement>& v1, const physics::measurements::measurement& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] - val);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator-(const physics::measurements::measurement& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(val - v1[i]);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator+=(const std::vector<physics::measurements::measurement>& v1, const std::vector<physics::measurements::measurement>& v2) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + v2[i]);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator+=(const std::vector<physics::measurements::measurement>& v1, const physics::measurements::measurement& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] + val);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator+=(const physics::measurements::measurement& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(val + v1[i]);
        return vec;
    }

    std::vector<physics::measurements::measurement> operator-=(const std::vector<physics::measurements::measurement>& v1, const std::vector<physics::measurements::measurement>& v2) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] - v2[i]);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator-=(const std::vector<physics::measurements::measurement>& v1, const physics::measurements::measurement& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] - val);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator-=(const physics::measurements::measurement& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(val - v1[i]);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator*(const std::vector<physics::measurements::measurement>& v1, const std::vector<physics::measurements::measurement>& v2) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] * v2[i]);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator*(const std::vector<physics::measurements::measurement>& v1, const physics::measurements::measurement& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] * val); 
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator*(const physics::measurements::measurement& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] * val); 
        return vec;
    }

    std::vector<physics::measurements::measurement> operator*(const std::vector<physics::measurements::measurement>& v1, const double& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] * val); 
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator*(const double& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] * val); 
        return vec;
    }

    std::vector<physics::measurements::measurement> operator/(const std::vector<physics::measurements::measurement>& v1, const std::vector<physics::measurements::measurement>& v2) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] / v2[i]);
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator/(const std::vector<physics::measurements::measurement>& v1, const physics::measurements::measurement& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] / val); 
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator/(const physics::measurements::measurement& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(val / v1[i]); 
        return vec;
    }

    std::vector<physics::measurements::measurement> operator/(const std::vector<physics::measurements::measurement>& v1, const double& val) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i] / val); 
        return vec;
    }
    
    std::vector<physics::measurements::measurement> operator/(const double& val, const std::vector<physics::measurements::measurement>& v1) {
        std::vector<physics::measurements::measurement> vec;
        vec.reserve(v1.size()); 
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(val / v1[i]); 
        return vec;
    }

    template <typename T>            
    std::vector<T> pow(const std::vector<T>& v1, const int& power) {
        std::vector<T> vec;
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i].pow(power));
        return vec;
    }

    template <typename T>            
    std::vector<T> square(const std::vector<T>& v1) {
        std::vector<T> vec;
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i].square());
        return vec;
    }
    
    template <typename T>            
    std::vector<T> cube(const std::vector<T>& v1) {
        std::vector<T> vec;
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i].cube());
        return vec;
    }
    
    template <typename T>            
    std::vector<T> root(const std::vector<T>& v1, const int& power) {
        std::vector<T> vec;
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i].root(power));
        return vec;
    }
    
    template <typename T>            
    std::vector<T> sqrt(const std::vector<T>& v1) {
        std::vector<T> vec;
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i].sqrt());
        return vec;
    }
    
    template <typename T>            
    std::vector<T> cbrt(const std::vector<T>& v1) {
        std::vector<T> vec;
        for (size_t i{}; i < v1.size(); i++) vec.emplace_back(v1[i].cbrt());
        return vec;
    }


    // operations amongs std::pair<std::vector<physics::measurements::measurement>, std::vector<physics::measurements::measurement>>
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator+(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v2) {
        return std::make_pair(v1.first + v2.first, v1.second + v2.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator-(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v2) {
        return std::make_pair(v1.first - v2.first, v1.second - v2.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator*(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v2) {
        return std::make_pair(v1.first * v2.first, v1.second * v2.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator/(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v2) {
        return std::make_pair(v1.first / v2.first, v1.second / v2.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator*(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const double& val) {
        return std::make_pair(v1.first * val, v1.second * val);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator*(const double& val, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(v1.first * val, v1.second * val);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator/(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const double& val) {
        return std::make_pair(v1.first / val, v1.second / val);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator/(const double& val, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(val / v1.first, val / v1.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator+(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const physics::measurements::measurement& meas) {
        return std::make_pair(v1.first + meas, v1.second + meas);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator+(const physics::measurements::measurement& meas, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(meas + v1.first, meas + v1.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator-(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const physics::measurements::measurement& meas) {
        return std::make_pair(v1.first - meas, v1.second - meas);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator-(const physics::measurements::measurement& meas, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(meas - v1.first, meas - v1.second);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator*(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const physics::measurements::measurement& meas) {
        return std::make_pair(v1.first * meas, v1.second * meas);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator*(const physics::measurements::measurement& meas, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(v1.first * meas, v1.second * meas);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator/(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                                const physics::measurements::measurement& meas) {
        return std::make_pair(v1.first / meas, v1.second / meas);
    }

    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> operator/(const physics::measurements::measurement& meas, 
                                                                                const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(meas / v1.first, meas / v1.second);
    }

    
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> pow(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                        const int& power) {
        return std::make_pair(pow(v1.first, power), pow(v1.second, power));
    }
    
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> square(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                            std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(square(v1.first), square(v1.second));
    }
    
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> cube(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(cube(v1.first), cube(v1.second));
    }
    
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> root(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1, 
                                                                        const int& power) {
        return std::make_pair(root(v1.first, power), root(v1.second, power));
    }
    
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> sqrt(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(sqrt(v1.first), sqrt(v1.second));
    }
    
    inline std::pair<std::vector<physics::measurements::measurement>, 
                    std::vector<physics::measurements::measurement>> cbrt(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                        std::vector<physics::measurements::measurement>>& v1) {
        return std::make_pair(cbrt(v1.first), cbrt(v1.second));
    }


} // namespace std


namespace physics {


    namespace tools {

        using namespace measurements; 


        // class expressing the coordinate of an object in a 3D system
        class position {

            private: 
                
                // =============================================
                // class members
                // =============================================
                
                std::vector<measurement> pos_ = std::vector<measurement>(3); 


            public: 

                // =============================================
                // constructors & destructor
                // =============================================
                
                explicit position(const measurement& x = (0. * units::SI::m), 
                                  const measurement& y = (0. * units::SI::m), 
                                  const measurement& z = (0. * units::SI::m)) noexcept {
                    pos_.emplace_back(x); 
                    pos_.emplace_back(y); 
                    pos_.emplace_back(z); 
                }

                position(const std::vector<measurement>& pos) noexcept : pos_{pos} {}

                position(const position& pos) noexcept : position(pos.get_position()) {}

                ~position() = default; 


                // =============================================
                // operators
                // =============================================

                inline position operator+(const position& other) const { return { pos_ + other.pos_ }; } 

                inline position operator-(const position& other) const { return { pos_ - other.pos_ }; } 

                inline position operator*(const position& other) const { return { pos_ * other.pos_ }; } 

                inline position operator*(const double& val) const { return { pos_ * val }; }

                inline position operator*(const int& val) const { return { pos_ * val }; }

                inline position operator/(const position& other) const { return { pos_ / other.pos_ }; } 

                inline position operator/(const double& val) const { return { pos_ / val }; }

                inline position operator/(const int& val) const { return { pos_ / val }; }


                // =============================================
                // set, get and print methods
                // =============================================

                inline void set_position(const measurement& x, const measurement& y, const measurement& z) { pos_ = {x, y, z}; }

                inline void set_position(const std::vector<measurement>& pos) { pos_ = pos; }
                
                inline void set_position(const position& pos) { pos_ = pos.get_position(); }

                inline void x(const measurement& x) { pos_[0] = x; }
                
                inline void y(const measurement& y) { pos_[1] = y; }
                
                inline void z(const measurement& z) { pos_[2] = z; }

                inline position as_position() const { return *this; }

                inline std::vector<measurement> get_position() const { return pos_; }

                inline measurement x() const { return pos_[0]; }
                
                inline measurement y() const { return pos_[1]; }
                
                inline measurement z() const { return pos_[2]; }             

                inline measurement magnitude() const { return (pos_[0].square() + pos_[1].square() + pos_[2].square()).sqrt(); }

                inline measurement distance(const position& other) {        
                    return ((pos_[0] - other.pos_[0]).square() + (pos_[1] - other.pos_[1]).square() + (pos_[2] - other.pos_[2]).square()).sqrt();
                }       

                inline measurement rho() const { return (pos_[0].square() + pos_[1].square()).sqrt(); }

                inline measurement phi() const { return measurement(atan2(pos_[1].value(), pos_[0].value()), units::SI::rad); }     

                inline measurement phi(const position& other) const { 
                    return measurement(atan2(other.pos_[1].value() - pos_[1].value(), other.pos_[0].value() - pos_[0].value()), units::SI::rad); 
                }

                inline measurement theta() const { 
                    if (pos_[2] == 0) return measurement(0, units::SI::rad);
                    return measurement(acos(pos_[2].value() / magnitude().value()), units::SI::rad); 
                }

                inline measurement theta(const position& other) {
                    if (other.pos_[2] == pos_[2]) return measurement(0, units::SI::rad);
                    return measurement(acos((other.pos_[2].value() - pos_[2].value()) / distance(other).value()), units::SI::rad); 
                }

                inline std::vector<double> direction() const {
                    return {cos(phi().value()), sin(phi().value()), pos_[2].value() / magnitude().value()};
                } 

                inline std::vector<double> direction(position& other) {
                    return {cos(phi(other).value()), sin(phi(other).value()), (other.pos_[2] - pos_[2]).value() / distance(other).value()};
                } 

                void print_position() const {
                    std::cout << "position = {\n";
                    for (auto i : pos_) i.print_measurement(); 
                    std::cout << "}\n"; 
                }


        }; // class position

        
        // class expressing the velocity of an object in a 3D system
        class velocity {

            private: 
                
                // =============================================
                // class members
                // =============================================
                
                std::vector<measurement> vel_; 
            

            public:  

                // =============================================
                // constructors
                // =============================================

                explicit velocity(const measurement& x = measurement(0., units::SI::mps), 
                                    const measurement& y = measurement(0., units::SI::mps), 
                                    const measurement& z = measurement(0., units::SI::mps)) noexcept {
                    vel_.reserve(3); 
                    vel_.emplace_back(x); 
                    vel_.emplace_back(y); 
                    vel_.emplace_back(z); 
                }

                explicit velocity(const double& x, const double& y, const double& z) noexcept {
                    vel_.reserve(3); 
                    vel_.emplace_back(x, units::SI::mps);
                    vel_.emplace_back(y, units::SI::mps);
                    vel_.emplace_back(z, units::SI::mps);
                }

                velocity(const std::vector<measurement>& pos) noexcept : vel_{pos} {}

                velocity(const velocity& pos) noexcept : velocity(pos.get_velocity()) {}

                ~velocity() = default; 


                // =============================================
                // operators
                // =============================================

                inline velocity operator+(const velocity& other) const { return velocity(vel_ + other.vel_); } 

                inline velocity operator-(const velocity& other) const { return velocity(vel_ - other.vel_); } 

                inline velocity operator*(const velocity& other) const { return velocity(vel_ * other.vel_); } 

                inline velocity operator*(const double& val) const { return velocity(vel_ * val); }
                
                inline velocity operator*(const int& val) const { return velocity(vel_ * val); }

                inline velocity operator/(const velocity& other) const { return velocity(vel_ / other.vel_); } 

                inline velocity operator/(const double& val) const { return velocity(vel_ / val); }

                inline velocity operator/(const int& val) const { return velocity(vel_ / val); }


                // =============================================
                // set & get & print methods
                // =============================================

                inline void set_velocity(const measurement& x, const measurement& y, const measurement& z) { 
                    vel_ = std::vector<measurement>{x, y, z}; 
                }

                inline void set_velocity(const double& x, const double& y, const double& z) { 
                    vel_ = std::vector<measurement>{measurement(x, units::SI::mps), 
                                                                    measurement(y, units::SI::mps), 
                                                                    measurement(z, units::SI::mps)}; 
                }

                inline void set_velocity(const std::vector<measurement>& pos) { vel_ = pos; }
                
                inline void set_velocity(const velocity& pos) { vel_ = pos.get_velocity(); }

                inline void x(const measurement& x) { vel_[0] = x; }
                
                inline void y(const measurement& y) { vel_[1] = y; }
                
                inline void z(const measurement& z) { vel_[2] = z; }

                inline void x(const double& x) { vel_[0] = measurement(x, vel_[0].units()); }

                inline void y(const double& y) { vel_[1] = measurement(y, vel_[1].units()); }

                inline void z(const double& z) { vel_[2] = measurement(z, vel_[2].units()); }

                inline velocity as_velocity() const { return *this; }

                inline std::vector<measurement> get_velocity() const { return vel_; }

                inline measurement get_vx() const { return vel_[0]; }
                
                inline measurement get_vy() const { return vel_[1]; }
                
                inline measurement get_vz() const { return vel_[2]; }             

                inline measurement magnitude() const { return (vel_[0].square() + vel_[1].square() + vel_[2].square()).sqrt(); }

                inline measurement phi() const { 
                    return measurement(atan2(vel_[1].value(), vel_[0].value()), units::SI::rad); 
                }     

                inline measurement theta() const { 
                    if (vel_[2] == 0) return measurement(0, units::SI::rad);
                    return measurement(acos(vel_[2].value() / magnitude().value()), units::SI::rad); 
                }

                inline std::vector<double> direction() const {
                    return {cos(phi().value()), sin(phi().value()), vel_[2].value() / magnitude().value()};
                } 

                void print_velocity() const {
                    std::cout << "velocity = {\n";
                    for (auto i : vel_) i.print_measurement(); 
                    std::cout << "}\n"; 
                }


        }; // class velocity
        

        inline position operator*(const double& val, position& pos) { return position(pos * val); }
        
        inline position operator*(const int& val, position& pos) { return position(pos * val); }
        
        inline velocity operator*(const double& val, velocity& vel) { return velocity(vel * val); }
        
        inline velocity operator*(const int& val, velocity& vel) { return velocity(vel * val); }


        // class expressing the position and velocity of a material point in a 3D system
        class material_point : public position, public velocity {

            public: 

                // =============================================
                // constructors and destructor
                // =============================================

                explicit material_point(const position& pos = position(), const velocity& vel = velocity()) noexcept : 
                    position(pos), velocity(vel) {}

                material_point(const std::pair<position, velocity>& pos_vel) noexcept : 
                    position(pos_vel.first), velocity(pos_vel.second) {}

                material_point(const material_point& other) noexcept : 
                    position(other.get_position()), velocity(other.get_velocity()) {}    
                
                ~material_point() = default;


                // =============================================
                // set, get and print methods
                // =============================================

                void set_pos_vel(const position& pos, const velocity& vel) {
                    set_position(pos); 
                    set_velocity(vel); 
                }

                void set_pos_vel(const std::pair<position, velocity>& pos_vel) {
                    set_position(pos_vel.first); 
                    set_velocity(pos_vel.second); 
                }               

                void set_pos_vel(const std::pair<std::vector<measurement>, std::vector<measurement>>& pos_vel) {
                    set_position(pos_vel.first); 
                    set_velocity(pos_vel.second); 
                }               

                inline std::pair<std::vector<measurement>, std::vector<measurement>> get_pos_vel() const { 
                    return std::make_pair(get_position(), get_velocity()); 
                }

                inline material_point as_material_point() const { return *this; }

                void print_pos_vel() const {
                    print_position(); 
                    print_velocity(); 
                }


        }; // class material point

        // class expressing the possible shapes for an object in a 3D system
        class shape {

            protected: 

                // =============================================
                // class member
                // =============================================

                measurement perimetre_, area_, volume_; 
                
            
            public: 

                // =============================================
                // virtual destructor
                // =============================================

                virtual ~shape() {}


                // =============================================
                // set and get methods
                // =============================================

                constexpr measurement perimetre() const { return perimetre_; }
                
                constexpr measurement area() const { return area_; }

                constexpr measurement volume() const { return volume_; }


                // =============================================
                // print methods
                // =============================================

                void print_perimetre() const { 
                    std::cout << "perimetre = "; 
                    perimetre_.print_measurement(); 
                }                        
                
                void print_area() const { 
                    std::cout << "area = "; 
                    area_.print_measurement(); 
                }        
                
                void print_volume() const { 
                    std::cout << "volume = "; 
                    volume_.print_measurement(); 
                }

        }; // class shape


        class sphere : public shape {

            protected:
                
                // =============================================
                // class member
                // =============================================

                measurement radius_; 


            public: 

                // =============================================
                // constructor and destructor
                // =============================================

                sphere(const measurement& radius) : radius_{radius} { 
                    perimetre_ = 2 * math::constants::pi * radius_;
                    area_ = 4 * math::constants::pi * radius_.square(); 
                    volume_ = 4. * math::constants::pi * radius_.cube() / 3.; 
                }
                
                ~sphere() {}


                // =============================================
                // set and get methods
                // =============================================

                constexpr void radius(const double& radius) { 
                    radius_ = radius; 
                    perimetre_ = 2 * math::constants::pi * radius_; 
                    area_ = 4 * math::constants::pi * radius_.square(); 
                    volume_ = 4. * math::constants::pi * radius_.cube() / 3.;  
                }

                constexpr measurement radius() const { return radius_; }

                void print_radius() const { 
                    std::cout << "radius = "; 
                    radius_.print_measurement(); 
                }

        }; // class sphere

        
    } // namespace tools 

    namespace objects {

    
        class mass : virtual public tools::material_point {

            protected: 

                // =============================================
                // class member
                // =============================================

                measurements::measurement mass_; 
                

            public: 

                // =============================================
                // constructors and destructor
                // =============================================
                
                explicit mass(const double& mass = 0., const units::unit& unit = units::SI::kg, const tools::position& pos = position(), const tools::velocity& vel = velocity()) noexcept : 
                    tools::material_point(pos, vel), mass_{mass, unit} {}

                explicit mass(const double& mass, const units::unit& unit, const std::pair<position, velocity>& pos_vel) noexcept :
                    tools::material_point(pos_vel), mass_{mass, unit} {}

                explicit mass(const measurements::measurement& mass, const tools::position& pos = position(), const tools::velocity& vel = velocity()) noexcept :
                    tools::material_point(pos, vel), mass_{mass} {}

                explicit mass(const measurements::measurement& mass, const std::pair<tools::position, tools::velocity>& pos_vel) noexcept :
                    tools::material_point(pos_vel), mass_{mass} {}

                mass(const mass& mass) noexcept : tools::material_point(mass.get_pos_vel()), mass_{mass.mass_} {}

                ~mass() = default; 


                // =============================================
                // operators
                // =============================================

                inline mass operator+(const mass& other) { return mass(mass_ + other.mass_, get_pos_vel()); }                    
                
                inline mass operator-(const mass& other) { return mass(mass_ - other.mass_, get_pos_vel()); }                    

                inline mass operator+=(const mass& other) { return mass(mass_ + other.mass_, get_pos_vel()); }                    
                
                inline mass operator-=(const mass& other) { return mass(mass_ - other.mass_, get_pos_vel()); }                    
                
                inline mass operator*(const mass& other) { return mass(mass_ * other.mass_, get_pos_vel()); }                    
                
                inline mass operator*(const double& val) { return mass(mass_ * val, get_pos_vel()); }  

                inline mass operator/(const mass& other) { return mass(mass_ / other.mass_, get_pos_vel()); }                    

                inline mass operator/(const double& val) { return mass(mass_ / val, get_pos_vel()); }  


                // =============================================
                // set, get and print methods
                // =============================================

                constexpr void set_mass(const double& mass) { mass_.value(mass); }
                
                constexpr measurements::measurement get_mass() const { return mass_; }

                inline mass as_mass() const { return *this; }

                void print_mass() const { 
                    std::cout << "mass = "; 
                    mass_.print_measurement(); 
                }
            

        }; // class mass


        class charge : virtual public tools::material_point {

            protected: 

                // =============================================
                // class member
                // =============================================

                measurements::measurement charge_; 
                

            public: 

                // =============================================
                // constructors and destructor
                // =============================================
                
                explicit charge(const double& charge = 0., const units::unit& unit = units::SI::derived::C, const tools::position& pos = position(), const tools::velocity& vel = velocity()) noexcept : 
                    material_point(pos, vel), charge_{charge, unit} {}

                explicit charge(const double& charge, const units::unit& unit, const std::pair<tools::position, tools::velocity>& pos_vel) noexcept :
                    material_point(pos_vel), charge_{charge, unit} {}

                explicit charge(const measurements::measurement& charge, const tools::position& pos = position(), const tools::velocity& vel = velocity()) noexcept :
                    material_point(pos, vel), charge_{charge} {}

                explicit charge(const measurements::measurement& charge, const std::pair<position, velocity>& pos_vel) noexcept :
                    material_point(pos_vel), charge_{charge} {}

                charge(const charge& charge) noexcept : material_point(charge.get_pos_vel()), charge_{charge.charge_} {}

                ~charge() = default; 


                // =============================================
                // operators
                // =============================================

                inline charge operator+(const charge& other) { return charge(charge_ + other.charge_, get_pos_vel()); }                    
                
                inline charge operator-(const charge& other) { return charge(charge_ - other.charge_, get_pos_vel()); }                    

                inline charge operator+=(const charge& other) { return charge(charge_ + other.charge_, get_pos_vel()); }                    
                
                inline charge operator-=(const charge& other) { return charge(charge_ - other.charge_, get_pos_vel()); }                    
                
                inline charge operator*(const charge& other) { return charge(charge_ * other.charge_, get_pos_vel()); }                    
                
                inline charge operator*(const double& val) { return charge(charge_ * val, get_pos_vel()); }  

                inline charge operator/(const charge& other) { return charge(charge_ / other.charge_, get_pos_vel()); }                    

                inline charge operator/(const double& val) { return charge(charge_ / val, get_pos_vel()); }  


                // =============================================
                // set, get and print methods
                // =============================================

                constexpr void set_charge(const double& charge) { charge_.value(charge); }
                
                constexpr measurements::measurement get_charge() const { return charge_; }

                inline charge as_charge() const { return *this; }

                void print_charge() const { 
                    std::cout << "charge = "; 
                    charge_.print_measurement(); 
                }
            

        }; // class charge
        

        class particle : public mass, public charge {

            public:  

                // =============================================
                // constructors & destructor
                // =============================================     

                explicit particle(const mass& mass = mass(), const charge& charge = charge(), const material_point& pos_vel = material_point()) noexcept : 
                    material_point(pos_vel), mass::mass(mass), charge::charge(charge) {}

                explicit particle(const mass& mass, const charge& charge, const tools::position& pos, const tools::velocity& vel) noexcept : 
                    material_point(pos, vel), mass::mass(mass), charge::charge(charge) {}

                explicit particle(const mass& mass, const charge& charge, const std::pair<position, velocity>& pos_vel) noexcept : 
                    material_point(pos_vel), mass::mass(mass), charge::charge(charge) {}

                explicit particle(const measurements::measurement& mass, const measurements::measurement& charge, const material_point& pos_vel = material_point()) noexcept : 
                    material_point(pos_vel), mass::mass(mass), charge::charge(charge) {}

                explicit particle(const measurements::measurement& mass, const measurements::measurement& charge, const tools::position& pos, const tools::velocity& vel) noexcept : 
                    material_point(pos, vel), mass::mass(mass), charge::charge(charge) {}

                explicit particle(const measurements::measurement& mass, const measurements::measurement& charge, const std::pair<position, velocity>& pos_vel) noexcept : 
                    material_point(pos_vel), mass::mass(mass), charge::charge(charge) {}

                particle(const particle& part) noexcept : particle(part.get_mass(), part.get_charge(), part.material_point::get_pos_vel()) {}

                ~particle() = default; 


                // =============================================
                // print methods
                // =============================================

                inline particle as_particle() const { return *this; }

                void print_particle() const {
                    std::cout << "particle : \n";
                    print_mass();
                    print_charge(); 
                    print_pos_vel();
                }


        }; // class particle


        class celestial_body : public mass {

            protected: 

                // =============================================
                // class members
                // =============================================

                const char* name_; 

                const char* celestial_type_; 


            public: 

                // =============================================
                // constructors and destructor
                // =============================================
                
                celestial_body(const char* type, const char* name, const mass& mass) : 
                    mass::mass(mass), name_{name}, celestial_type_{type} {}

                ~celestial_body() {}


                // =============================================
                // set and get methods
                // =============================================

                constexpr void name(const char* name) { name_ = name; }

                const char* name() const { return name_; }

                const char* celestial_type() const { return celestial_type_; }

                
                // =============================================
                // print methods
                // =============================================

                void print_name() const { std::cout << "name = " << name() << "\n"; }
                
                void print_celestial_type() const { std::cout << "type = " << celestial_type() << "\n"; }

                void print_celestial_body() const {
                    std::cout << "Celestial body :\n"; 
                    print_celestial_type(); 
                    print_name(); 
                    print_mass();
                    print_position();
                }

        }; // class celestial_body


        class planet : public celestial_body, public tools::sphere {

            protected:

                // =============================================
                // class members
                // =============================================
                
                position aphelion_pos_{}, perihelion_pos_{}; 
                velocity aphelion_vel_{}, perihelion_vel_{}; 
                measurements::measurement period_{};

            public: 

                // =============================================
                // constructors and destructor
                // =============================================

                planet(const char* name, const mass& mass, const measurements::measurement& radius) : celestial_body("planet", name, mass), tools::sphere(radius) {}

                planet(const char* name, const mass& mass, 
                       const measurements::measurement& radius, const measurements::measurement& period,
                       const position& aphelion_pos = position(), const position& perihelion_pos = position(), 
                       const tools::velocity& aphelion_vel = velocity(), const tools::velocity& perihelion_vel = velocity()) :
                    celestial_body("planet", name, mass), tools::sphere(radius),
                    aphelion_pos_{aphelion_pos}, perihelion_pos_{perihelion_pos}, 
                    aphelion_vel_{aphelion_vel}, perihelion_vel_{perihelion_vel}, 
                    period_{period} {}

                
                // =============================================
                // get methods
                // =============================================

                inline position pos_aphelion() { return aphelion_pos_; }

                inline position pos_perihelion() { return perihelion_pos_; }                
                
                inline velocity vel_aphelion() { return aphelion_vel_; }

                inline velocity vel_perihelion() { return perihelion_vel_; }

                inline measurements::measurement period() const { return period_; }


                // =============================================
                // print methods
                // =============================================

                void print_period() const { 
                    std::cout << "period = "; 
                    period_.print_measurement();
                }

                void print_planet() const {
                    print_celestial_body(); 
                    print_radius(); 
                    print_period(); 
                }


        }; // class planet
    

    } // namespace objects


    namespace solar_system {

        using objects::planet; 
        using objects::mass;
        using namespace measurements;

        planet Sun("Sun", mass(1.98844E30 * kg), 695700 * km, 0.0 * s); 
        
        planet Mercury("Mercury", mass(0.33010E24 * kg), 2440.5 * km, 87.969 * day,
                                    tools::position(69.818E6 * km, 0. * km, 0. * km), tools::position(-46E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 38.86 * kmps, 0. * kmps), tools::velocity(0. * kmps, -58.98 * kmps, 0. * kmps)); 

        planet Venus("Venus", mass(4.8673E24 * kg), 6051.8 * km, 224.701 * day,
                                    tools::position(108.941E6 * km, 0. * km, 0. * km), tools::position(-107.480E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 34.79 * kmps, 0. * kmps), tools::velocity(0. * kmps, -35.26 * kmps, 0. * kmps)); 

        planet Earth("Earth", mass(5.9722E24 * kg), 6378.137 * km, 365.256 * day,
                                    tools::position(152.100E6 * km, 0. * km, 0. * km), tools::position(-147.095E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 29.2911 * kmps, 0. * kmps), tools::velocity(0. * kmps, -30.2865 * kmps, 0. * kmps)); 
       
        planet Mars("Mars", mass(0.64169E24 * kg), 3396.2 * km, 686.980 * day,
                                    tools::position(249.261E6 * km, 0. * km, 0. * km), tools::position(-206.650E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 21.97 * kmps, 0. * kmps), tools::velocity(0. * kmps, -26.50 * kmps, 0. * kmps)); 

        planet Jupiter("Jupiter", mass(1898.13E24 * kg), 71492 * km, 4332.589 * day,
                                    tools::position(816.363E6 * km, 0. * km, 0. * km), tools::position(-740.595E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 12.44 * kmps, 0. * kmps), tools::velocity(0. * kmps, -13.72 * kmps, 0. * kmps)); 

        planet Saturn("Saturn", mass(568.32E24 * kg), 60268 * km, 10759.22 * day,
                                    tools::position(108.941E6 * km, 0. * km, 0. * km), tools::position(-107.480E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 9.09 * kmps, 0. * kmps), tools::velocity(0. * kmps, -10.18 * kmps, 0. * kmps)); 

        planet Uranus("Uranus", mass(86.811E24 * kg), 25559* km, 30685.4 * day,
                                    tools::position(3001.390E6 * km, 0. * km, 0. * km), tools::position(-2732.696E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 6.49 * kmps, 0. * kmps), tools::velocity(0. * kmps, -7.11 * kmps, 0. * kmps)); 
       
        planet Neptune("Neptune", mass(102.409E24 * kg), 24764 * km, 60189 * day,
                                    tools::position(4558.857E6 * km, 0. * km, 0. * km), tools::position(-4471.050E6 * km, 0. * km, 0. * km),
                                    tools::velocity(0 * kmps, 5.37 * kmps, 0. * kmps), tools::velocity(0. * kmps, -5.50 * kmps, 0. * kmps)); 


    } // namespace solar_system


    // namespace defining some usefull constants
    namespace constants {

        objects::mass e_mass = objects::mass(9.109383701528e-31, units::SI::kg);     
        objects::mass p_mass = objects::mass(1.672621923695e-27, units::SI::kg);
        objects::charge e_charge = objects::charge(-1.602176634e-19, units::SI::derived::C);
        objects::charge p_charge = objects::charge(1.602176634e-19, units::SI::derived::C); 

    } // namespace constants


    namespace objects {


        particle electon = particle(constants::e_mass, constants::e_charge); 

        particle proton = particle(constants::p_mass, constants::p_charge);


        class time {

            public:     

                // =============================================
                // class members
                // =============================================     

                measurements::measurement time_; 


                // =============================================
                // constructor and destructor
                // =============================================   

                explicit constexpr time(const double& time = 0.0, const units::unit& unit = units::SI::s) noexcept : time_(time, unit) {}
                
                ~time() = default; 

                
                // =============================================
                // set & get & print methods
                // =============================================   

                constexpr void set_time(const double& t) { time_.value(t); }

                constexpr void set_time(const measurements::measurement& t) { time_ = t; }
                                    
                constexpr measurements::measurement get_time() const { return time_; }

                constexpr time as_time() const { return *this; }

                constexpr void increase_time(const measurements::measurement& t) { time_ = time_ + t; } 

                constexpr void reset_time() { time_.value(0.0); }

                inline void print_time() const { time_.print_measurement(); }


        }; // class time


        // class for keeping track of the inexorable passage of time
        class timer {
            
            protected: 

                // =============================================
                // class members
                // =============================================     
                
                std::chrono::time_point<std::chrono::system_clock> start_, pause_;

                units::unit unit_; 

            
            public:

                // =============================================
                // constructor and destructor
                // =============================================   

                explicit constexpr timer(const units::unit& unit = units::SI::s) noexcept : unit_{unit} {}

                ~timer() = default;

                
                // =============================================
                // timer methods
                // =============================================   

                inline void start() { start_ = std::chrono::system_clock::now(); }

                inline void pause() { pause_ = std::chrono::system_clock::now(); }
                
                void print() const { 
                    std::cout << "elapsed time = " << static_cast<std::chrono::duration<double>>(pause_ - start_).count() << " "; 
                    unit_.print_unit(); 
                    std::cout << "\n";
                }


        }; // class timer


    } // namespace objects


} // namespace physics


namespace math {


    namespace tools {

        using namespace algebra; 


        class ode_solver : public physics::objects::time {

                protected: 

                    // =============================================
                    // class members
                    // =============================================

                    std::pair<std::vector<physics::measurements::measurement>, 
                                std::vector<physics::measurements::measurement>> m_df; 
                    

                public: 

                    // =============================================
                    // virtual destructor
                    // =============================================

                    virtual ~ode_solver() {}


                    // =============================================
                    // virtual diff methods
                    // =============================================

                    virtual std::pair<std::vector<physics::measurements::measurement>,
                                        std::vector<physics::measurements::measurement>> diff(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                            std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                            const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) = 0; 


                    // =============================================
                    // integration methods
                    // =============================================

                    inline std::pair<std::vector<physics::measurements::measurement>,
                                        std::vector<physics::measurements::measurement>> euler(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                            std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                            const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) {
                        increase_time(h); 
                        return pos_vel + h * diff(pos_vel, h);
                    } 

                    std::pair<std::vector<physics::measurements::measurement>,
                                std::vector<physics::measurements::measurement>> euler_modified(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                                std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                                const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) {
                        std::pair<std::vector<physics::measurements::measurement>,
                                    std::vector<physics::measurements::measurement>> appo = pos_vel + h * diff(pos_vel, h); 
                        increase_time(h); 
                        return pos_vel + h * (appo + diff(pos_vel + h * appo, h)) / 2.; 
                    }

                    std::pair<std::vector<physics::measurements::measurement>,
                                std::vector<physics::measurements::measurement>> rk4(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                    std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                    const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) {
                        std::pair<std::vector<physics::measurements::measurement>,
                                    std::vector<physics::measurements::measurement>> k1{}, k2{}, k3{}, k4{}; 
                        k1 = diff(pos_vel, h); 
                        k2 = diff(pos_vel + k1 * h / 2., h);
                        k3 = diff(pos_vel + k2 * h / 2., h);
                        k4 = diff(pos_vel + k3 * h / 2., h);      
                        increase_time(h); 
                        return (pos_vel + (k1 + k2 * 2. + k3 * 2. + k4) * (h / 6.)); 
                    } 


            }; // class ode_solver


    } // namespace tools


} // namespace math


namespace physics {


    namespace objects {


        namespace oscillators {


            class harmonic : public math::tools::ode_solver {

                private: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::measurement omega_;


                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    harmonic(const measurements::measurement& omega, const measurements::measurement& time = measurements::measurement(0.0, units::SI::s)) : omega_{omega} { time::time_ = time; }
                    
                    ~harmonic() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    constexpr void set_omega(const measurements::measurement& omega) { omega_ = omega; }

                    constexpr measurements::measurement get_omega() const { return omega_; }

                    void print_omega() const {
                        std::cout << "omega = "; 
                        omega_.print_measurement();
                    }


                    // =============================================
                    // diff methods
                    // =============================================

                    inline std::pair<std::vector<measurements::measurement>, 
                                        std::vector<measurements::measurement>> diff(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                    std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                    const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) override {
                        return std::make_pair(pos_vel.second, - omega_.square() * pos_vel.first);
                    }


            }; // class harmonic


            class forced : public math::tools::ode_solver {

                private: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::measurement omega_0_, omega_1_;


                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    forced(const measurements::measurement& omega0, 
                           const measurements::measurement& omega1, 
                           const measurements::measurement& time = measurements::measurement(0.0, units::SI::s)) : 
                        omega_0_{omega0}, omega_1_{omega1} { time::time_ = time; }
                    
                    ~forced() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    constexpr void set_omega0(const measurements::measurement& omega) { omega_0_ = omega; }

                    constexpr void set_omega1(const measurements::measurement& omega) { omega_1_ = omega; }

                    constexpr measurements::measurement get_omega0() const { return omega_0_; }

                    constexpr measurements::measurement get_omega1() const { return omega_1_; }

                    void print_omega0() const {
                        std::cout << "omega0 = "; 
                        omega_0_.print_measurement();
                    }

                    void print_omega1() const {
                        std::cout << "omega1 = "; 
                        omega_1_.print_measurement();
                    }


                    // =============================================
                    // diff methods
                    // =============================================

                    inline std::pair<std::vector<measurements::measurement>, 
                                        std::vector<measurements::measurement>> diff(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                    std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                    const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) override {
                        return std::make_pair(pos_vel.second, - omega_0_.square() * pos_vel.first + std::sin(omega_1_.value() * time::time_.value()));
                    }


            }; // class forced


            class damped : public math::tools::ode_solver {

                private: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::measurement omega_, alpha_;


                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    damped(const measurements::measurement& omega, 
                           const measurements::measurement& alpha, 
                           const measurements::measurement& time = measurements::measurement(0.0, units::SI::s)) : 
                        omega_{omega}, alpha_{alpha} { time::time_ = time; }
                    
                    ~damped() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    constexpr void set_omega(const measurements::measurement& omega) { omega_ = omega; }

                    constexpr void set_alpha(const measurements::measurement& alpha) { alpha_ = alpha; }

                    constexpr measurements::measurement get_omega() const { return omega_; }

                    constexpr measurements::measurement get_alpha() const { return alpha_; }

                    void print_omega() const {
                        std::cout << "omega = "; 
                        omega_.print_measurement();
                    }

                    void print_alpha() const {
                        std::cout << "alpha = "; 
                        alpha_.print_measurement();
                    }


                    // =============================================
                    // diff methods
                    // =============================================

                    inline std::pair<std::vector<measurements::measurement>, 
                                        std::vector<measurements::measurement>> diff(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                    std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                    const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) override {
                        return std::make_pair(pos_vel.second, - omega_.square() * pos_vel.first - alpha_ * pos_vel.second);
                    }


            }; // class damped


            class forced_damped : public math::tools::ode_solver {

                private: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::measurement omega_0_, omega_1_, alpha_;


                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    forced_damped(const measurements::measurement& omega0, 
                                  const measurements::measurement& omega1, 
                                  const measurements::measurement& alpha,
                                  const measurements::measurement& time = measurements::measurement(0.0, units::SI::s)) : 
                        omega_0_{omega0}, omega_1_{omega1}, alpha_{alpha} { time::time_ = time; }
                    
                    ~forced_damped() {}


                    // =============================================
                    // set and get methods
                    // =============================================

                    constexpr void set_omega0(const measurements::measurement& omega) { omega_0_ = omega; }

                    constexpr void set_omega1(const measurements::measurement& omega) { omega_1_ = omega; }

                    constexpr void set_alpha(const measurements::measurement& alpha) { alpha_ = alpha; }

                    constexpr measurements::measurement get_omega0() const { return omega_0_; }

                    constexpr measurements::measurement get_omega1() const { return omega_1_; }

                    constexpr measurements::measurement get_alpha() const { return alpha_; }

                    void print_omega0() const {
                        std::cout << "omega0 = "; 
                        omega_0_.print_measurement();
                    }

                    void print_omega1() const {
                        std::cout << "omega1 = "; 
                        omega_1_.print_measurement();
                    }

                    void print_alpha() const {
                        std::cout << "alpha = "; 
                        alpha_.print_measurement();
                    }


                    // =============================================
                    // diff methods
                    // =============================================

                    inline std::pair<std::vector<measurements::measurement>, 
                                        std::vector<measurements::measurement>> diff(const std::pair<std::vector<physics::measurements::measurement>, 
                                                                                                    std::vector<physics::measurements::measurement>>& pos_vel, 
                                                                                    const physics::measurements::measurement& h = physics::measurements::measurement(0.001, physics::units::SI::s)) override {
                        return std::make_pair(pos_vel.second, - omega_0_.square() * pos_vel.first + std::sin(omega_1_.value() * time::time_.value()) - alpha_ * pos_vel.second);
                    }


            }; // class forced_damped


        } // namespace oscillators


    } // namespace objects


} // namespace physics

