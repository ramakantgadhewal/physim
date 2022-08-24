

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physim(namespace) containing the basic tools for computational physics. 
// last updated:    24/08/2022


#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <fstream>
#include <iomanip> 
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>


namespace physim {


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
            constexpr T square(const T& value) { return value * value; }

            // generate the cubic power of a value
            template<typename T> 
            constexpr T cube(const T& value) { return value * value * value; }

            // generate small integer powers of a value (1, 0, -1)
            template<typename T> 
            constexpr T pow_small(const T& value, const int& power) { 
                return (power == 1) ? value : ((power == -1) ? T(1.0) / value : T(1.0));
            }

            // generate _root power of a value
            template<typename T> T root(const T& value, const int& power) {
                switch (power) {
                    case 0: 
                        return T{1.0};
                    case 1: 
                        return value;
                    case -1: 
                        return T{1.0} / value;
                    case 2: 
                        if (value < T{}) { return constants::invalid_conversion; }
                        else return std::sqrt(value);
                    case -2:
                        if (value < T{}) { return constants::invalid_conversion; }
                        else return std::sqrt(T{1.0} / value);
                    case 3: 
                        return std::cbrt(value);
                    case -3: 
                        return std::cbrt(T{1.0} / value);
                    case 4: 
                        if (value < T{}) { return constants::invalid_conversion; }
                        else return std::sqrt(std::sqrt(value));
                    case -4: 
                        if (value < T{}) { return constants::invalid_conversion; }
                        else return std::sqrt(std::sqrt(T{1.0} / value));
                    default:
                        if (value < T{} && power % 2 == 0) { return constants::invalid_conversion; }
                        else return std::pow(value, T{1.0} / static_cast<T>(power));
                }
            }
            
            // generate the square root power of a value
            template <typename T> 
            T sqrt(const T& value) { return root(value, 2); }

            // generate the cubic root power of a value
            template <typename T> 
            T cbrt(const T& value) { return root(value, 3); }


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
                    result = 1.0 / static_cast<double>(i+1) * (static_cast<double>(i) * result + math::algebra::square(v[i]) + static_cast<double>(i) * math::algebra::square(old_average)) - math::algebra::square(average);
                }
                return result;
            }

            // standard deviation of an std::vector
            template <typename K>
            double sd(const std::vector<K>& v) {
                return math::algebra::sqrt(variance(v));
            }

            // standard deviation of mean of an std::vector
            template <typename K>
            double sdom(const std::vector<K>& v) {
                return sd(v) / math::algebra::sqrt(v.size());
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
                
                    virtual double eval(const double& x) const = 0; 
                
                    int signum(const double& x) const { return (eval(x) == 0. ? 0. : (eval(x) > 0 ? 1. : -1)); }


                    // =============================================
                    // print methods
                    // =============================================
                
                    void print_eval(const double& x, const double& precision = 1.e-6) const {
                        std::cout << "f (" << x << ") = " << std::setprecision(precision) << eval(x) << "\n"; 
                    }
                
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
                
                    functor(const char& op, function_base* f, function_base* g) {
                        op_ = op;
                        f_ = f; 
                        g_ = g;
                    }
                    
                    ~functor() {}
                
                
                    // =============================================
                    // eval methods
                    // =============================================
                
                    double eval(const double& x) const override {
                        switch (op_) {
                            case '+':
                                return f_->eval(x) + g_->eval(x);
                            case '-': 
                                return f_->eval(x) - g_->eval(x);
                            case '*':
                                return f_->eval(x) * g_->eval(x);
                            case '/': 
                                return f_->eval(x) / g_->eval(x);
                            case '^':
                                return std::pow(f_->eval(x), g_->eval(x));
                            case 'c':
                                return f_->eval(g_->eval(x));
                            default:
                                break;
                        }
                        return NAN;
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
            
                    explicit line(const double& m = 1, const double& q = 0) noexcept : m_{m}, q_{q} {}
                    
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
            
                    explicit quadratic(const double& a, const double& b, const double& c) noexcept :
                        a_{a}, b_{b}, c_{c}, delta_{math::algebra::square(b_) - 4 * a_ * c_} {}
                    
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

                    std::pair<double, double> roots() const {
                        if (delta_ == 0) return std::make_pair(- b_ / (2 * a_), - b_ / (2 * a_));
                        else if (delta_ > 0) return std::make_pair(- b_ - algebra::sqrt(delta_) / (2 * a_), (- b_ + algebra::sqrt(delta_)) / (2 * a_));
                        else std::cout << "\nThere are not real solutions...\n\n"; // note: create a complex class then come back here
                        return std::make_pair(NAN, NAN); 
                    }


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
                
                    explicit cubic(const double& a, const double& b, const double& c, const double& d) noexcept :
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
                        return a_ * math::algebra::cube(x) + b_ * math::algebra::square(x) + c_ * x + d_; 
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
                
                    explicit square_root(const double& c = 1) noexcept : c_{c} {}
                
                    ~square_root() {}
                
                
                    // =============================================
                    // set & get methods
                    // =============================================
                
                    constexpr void c(const double& c) { c_ = c; }

                    constexpr double c() const { return c_; }

                
                    // =============================================
                    // eval methods
                    // =============================================

                    double eval(const double& x) const override { return c_ * math::algebra::sqrt(x); }
                            
            
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
                
                    explicit cubic_root(double c = 1) noexcept : c_{c} {}
                
                    ~cubic_root() {} 
                
                
                    // =============================================
                    // set & get methods
                    // =============================================
                
                    constexpr void c(double c) { c_ = c; }

                    constexpr double c() const { return c_; }

                
                    // =============================================
                    // eval methods
                    // =============================================

                    double eval(const double& x) const override { return c_ * math::algebra::cbrt(x); }
                            
            
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
            
                    explicit exponential(const double& base = math::constants::e, const double& c1 = 1, const double& c2 = 1) noexcept :
                        base_{base}, c1_{c1}, c2_{c2} {}

                    explicit exponential(const double& c1 = 1, const double& c2 = 1) noexcept :
                        base_{math::constants::e}, c1_{c1}, c2_{c2} {}

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
                        if (base_ != math::constants::e) std::cout << base_ << "^"; 
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
            
                    explicit logarithm(const double& base = math::constants::e, const double& c1 = 1, const double& c2 = 1) noexcept : 
                        base_{base}, c1_{c1}, c2_{c2} {}

                    explicit logarithm(const double& c1 = 1, const double& c2 = 1) noexcept : 
                        base_{math::constants::e}, c1_{c1}, c2_{c2} {}

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
                        if (base_ != math::constants::e) std::cout << "_( " << base_ << ")";
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
                    
                    // y = c1 * sin (c2 * x)
                    double c1_, c2_; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
                    explicit sine(const double& c1 = 1, const double& c2 = 1) noexcept : c1_{c1}, c2_{c2} {} 
                
                    ~sine() {} 

                
                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void c1(const double& c1) { c1_ = c1; }
                
                    constexpr void c2(const double& c2) { c2_ = c2; }
                
                    constexpr double c1() const { return c1_; }

                    constexpr double c2() const { return c2_; }
                
                
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
                    
                    // y = c1 * cos (c2 * x)
                    double c1_, c2_; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
            
                    explicit cosine(const double& c1 = 1, const double& c2 = 1) noexcept : c1_{c1}, c2_{c2} {} 
                
                    ~cosine() {} 

                
                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void c1(const double& c1) { c1_ = c1; }
                
                    constexpr void c2(const double& c2) { c2_ = c2; }
                
                    constexpr double c1() const { return c1_; }

                    constexpr double c2() const { return c2_; }
                

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
                    
                    // y = c1 * tan (c2 * x)
                    double c1_, c2_; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
                    explicit tangent(const double& c1 = 1, const double& c2 = 1) noexcept : c1_{c1}, c2_{c2} {} 
                
                    ~tangent() {} 

                
                    // =============================================
                    // set methods
                    // =============================================
            
                    constexpr void c1(const double& c1) { c1_ = c1; }
                
                    constexpr void c2(const double& c2) { c2_ = c2; }
                
                
                    // =============================================
                    // get methods
                    // =============================================

                    constexpr double c1() const { return c1_; }

                    constexpr double c2() const { return c2_; }
                

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
                    
                    unsigned int a_, c_, m_, seed_;


                public: 
            
                    // =============================================
                    // constructors & destructor
                    // =============================================
                    
                    random_generator() :  a_{1664525}, c_{1013904223}, m_{static_cast<unsigned int>(std::pow(2, 31))} {}

                    random_generator(const unsigned int& seed) : a_{1664525}, c_{1013904223}, m_{static_cast<unsigned int>(std::pow(2, 31))}, seed_{seed} {}

                    ~random_generator() {}

                    // =============================================
                    // set & get methods
                    // =============================================

                    constexpr void a(const unsigned int& a) { a_ = a; }
                    
                    constexpr void c(const unsigned int& c) { c_ = c; }
                    
                    constexpr void m(const unsigned int& m) { m_ = m; }
                    
                    constexpr void seed(const unsigned int& seed) { seed_ = seed; }

                    constexpr unsigned int a() const { return a_; }
                    
                    constexpr unsigned int c() const { return c_; }
                    
                    constexpr unsigned int m() const { return m_; }
                    
                    constexpr unsigned int seed() const { return seed_; }


                    // =============================================
                    // distributions methods
                    // =============================================

                    constexpr double rand(const double& min = 0., const double& max = 1.) {
                        seed(static_cast<unsigned int>((a_ * seed_ + c_) % m_)); 
                        return min + std::fabs(max - min) * seed_ / m_; 
                    }

                    constexpr double exp(const double& mean) {
                        return - std::log(1 - rand()) / mean; 
                    }

                    double gauss_box_muller(const double& mean, const double& sigma) {
                        double s{rand()}, t{rand()}; 
                        double x = math::algebra::sqrt(-2 * std::log(s)) * std::cos(2 * math::constants::pi * t);
                        return mean + sigma * x;
                    }

                    double gauss_accept_reject(const double& mean, const double& sigma) {
                        double x{}, y{}, g{}; 
                        while (true) {
                            x = rand(-5., 5.); 
                            y = rand(); 
                            g = exp(math::algebra::square(x) / 2); 
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

                    integral() {}
                    
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
                        print_value(); 
                        print_error(); 
                    }
                    
                    
                    // =============================================
                    // integration methods
                    // =============================================

                    void midpoint(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                        begin_integration(a, b, n); 
                        for (unsigned int i{}; i < steps_; i++) { sum_ += (f.eval(a_ + (i + 0.5) * h_)); }
                        integral_ = sum_ * h_; 
                    }

                    void midpoint_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                        begin_integration(a, b, 1); 
                        while (true) {
                            old_integral_ = integral_; 
                            midpoint(a_, b_, f, steps_ * 2);
                            error_ = 4 * std::fabs(integral_ - old_integral_) / 3;
                            if (error_ < prec) break;
                        }    
                    }
                    
                    void trapexoid(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                        begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 2.);
                        for (unsigned int i{1}; i < steps_; i++) sum_ += f.eval(a_ + i * h_); 
                        integral_ = sum_ * h_; 
                    } 

                    void trapexoid_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                        begin_integration(a, b, 2, f.eval(a) + f.eval(b) / 2. + f.eval((a + b) / 2.)); 
                        while (true) {
                            old_integral_ = integral_; 
                            trapexoid(a_, b_, f, steps_ * 2);
                            error_ = 4 * std::fabs(integral_ - old_integral_) / 3;
                            if (error_ < prec) break;
                        }
                    }

                    void simpson(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                        if (n % 2 == 0) begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 3.);
                        else begin_integration(a, b, n + 1);  
                        for (unsigned int i{1}; i < steps_; i++) sum_ += 2 * (1 + i % 2) * (f.eval(a_ + i * h_)) / 3.; 
                        integral_ = sum_ * h_; 
                    }

                    void simpson_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                        begin_integration(a, b, 2, (f.eval(a) + f.eval(b)) / 3.); 
                        while (true) {
                            old_integral_ = integral_; 
                            simpson(a_, b_, f, steps_ * 2);
                            error_ = 16 * std::fabs(integral_ - old_integral_) / 15;
                            if (error_ < prec) break; 
                        }
                    }

                    void mean(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                        begin_integration(a, b, n); 
                        for (unsigned int i{}; i < n; i ++) sum_ += f.eval(rg_.rand(a, b)); 
                        integral_ = (b_ - a_) * sum_ / steps_; 
                    }

                    void mean_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                        std::vector<double> k{};
                        for (unsigned i{}; i < 10000; i++) {
                            mean(a, b, f);
                            k.emplace_back(integral_); 
                        }
                        double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                        mean(a, b, f, static_cast<unsigned int>(math::algebra::square(k_mean / prec))); 
                    }
            
                    void hit_or_miss(const double& a, const double& b, const functions::function_base& f, const double& fmax, unsigned int n = 1000) {
                        begin_integration(a, b, n); 
                        double x{}, y{}; 
                        unsigned int hits{};
                        for (unsigned int i{}; i < n; i ++) {
                            x = rg_.rand(a, b); 
                            y = rg_.rand(0., fmax);  
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
                        unsigned int N = static_cast<unsigned int>(math::algebra::square(k_mean / prec)); 
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

                    // check if the units have the same base unit 
                    constexpr bool has_same_base(const unit_data& other) const { return *this == other; }


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
                            (second_ * power) + root_Hertz_modifier(power),
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
                                                                    second_ / power,
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
                    
                    constexpr int metre() const { return metre_; }

                    constexpr int second() const { return second_; }

                    constexpr int kg() const { return kilogram_; }
                                        
                    constexpr int ampere() const { return ampere_; }
                    
                    constexpr int kelvin() const { return kelvin_; }
                    
                    constexpr int mole() const { return mole_; }
                    
                    constexpr int candela() const { return candela_; }

                    constexpr unit_data data() const { return *this; }

                    constexpr void print() const {
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
                    constexpr int root_Hertz_modifier(const int& power) const {
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
                    
                    constexpr unit_prefix() noexcept : multiplier_{1.0} {}
                    
                    constexpr unit_prefix(const double& mult) noexcept : multiplier_{mult} {} 


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
                    inline unit_prefix root(const int& power) const { return { math::algebra::root(multiplier_, power) }; }

                    inline unit_prefix sqrt() const { return { math::algebra::root(multiplier_, 2) }; }
                    
                    inline unit_prefix cbrt() const { return { math::algebra::root(multiplier_, 3) }; }
                    
                    
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
                    constexpr void print() const { 
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
                    // constructors
                    // ============================================= 
                    
                    // default constructor
                    constexpr unit() noexcept : unit_data(), unit_prefix() {}

                    // constructor from unit_data 
                    constexpr unit(const unit_data& data) noexcept : unit_data(data) {}

                    // constructor from unit_data and unit_prefix
                    constexpr unit(const unit_data& data, const unit_prefix& prefix) noexcept : unit_data(data), unit_prefix(prefix) {}

                    // constructor from unit_prefix and unit_data
                    constexpr unit(const unit_prefix& prefix, const unit_data& data) noexcept : unit_data(data), unit_prefix(prefix) {}

                    // constructor from unit_data and a double
                    constexpr unit(const unit_data& data, const double& mult) noexcept : unit_data(data), unit_prefix(mult) {}

                    // constructor from double with a unit_data
                    constexpr unit(const double& mult, const unit_data& data) noexcept : unit_data(data), unit_prefix(mult) {}

                    // constructor from unit
                    constexpr unit(const unit& other) noexcept : unit_data(other.data()), unit_prefix(other.prefix()) {} 

                    // constructor from unit and unit_prefix
                    constexpr unit(const unit& other, const unit_prefix& prefix) noexcept : unit_data(other.data()), unit_prefix(prefix) {} 

                    // constructor from unit_prefix and unit
                    constexpr unit(const unit_prefix& prefix, const unit& other) noexcept : unit_data(other.data()), unit_prefix(prefix) {} 

                    // constructor from unit and a double
                    constexpr unit(const unit& other, const double& multiplier) noexcept : unit_data(other.data()), unit_prefix(multiplier) {} 

                    // constructor from double with a unit
                    constexpr unit(const double& multiplier, const unit& other) noexcept : unit_data(other.data()), unit_prefix(multiplier) {} 


                    // =============================================
                    // operators
                    // ============================================= 

                    // multiply with another unit
                    constexpr unit operator*(const unit& other) const { 
                        return { unit_data::operator*(other.data()), unit_prefix::operator*(other.prefix()) }; 
                    }

                    // division operator
                    constexpr unit operator/(const unit& other) const { 
                        return { unit_data::operator/(other.data()), unit_prefix::operator/(other.prefix()) }; 
                    }                    

                    // equality operator
                    constexpr bool operator==(const unit& other) const {
                        if (unit_data::operator!=(other.data())) { return false; }
                        if (unit_prefix::operator==(other.prefix())) { return true; }
                        return math::tools::compare_round_equals(multiplier(), other.multiplier());
                    }

                    // equality operator
                    constexpr bool operator!=(const unit& other) const { return !operator==(other); }


                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a unit
                    constexpr unit inv() const { return { unit_data::inv(), unit_prefix::inv() }; }

                    // take a unit to a power
                    constexpr unit pow(const int& power) const { return { unit_data::pow(power), unit_prefix::pow(power) }; }

                    constexpr unit square() const { return { unit_data::square(), unit_prefix::square() }; }

                    constexpr unit cube() const { return { unit_data::cube(), unit_prefix::cube() }; }

                    // take a unit to a root power
                    inline unit root(const int& power) const { return { unit_data::root(power), unit_prefix::root(power) }; }

                    inline unit sqrt() const { return { unit_data::sqrt(), unit_prefix::sqrt() }; }

                    inline unit cbrt() const { return { unit_data::cbrt(), unit_prefix::cbrt() }; }


                    // =============================================
                    // checks
                    // ============================================= 

                    // test for exact numerical equivalence
                    constexpr bool is_exactly_the_same(const unit& other) const { 
                        return unit_data::operator==(other.data()) && unit_prefix::operator==(other.prefix());
                    }


                    // =============================================
                    // get methods
                    // ============================================= 

                    // get the unit
                    constexpr unit as_unit() const { return *this; }

                    // print the unit
                    constexpr void print() const {
                        unit_prefix::print(); 
                        unit_data::print();
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
            

            // derived_units
            namespace SI_derived {

                unit hertz(unit_data(0, -1, 0, 0, 0, 0, 0));
                unit Hz = hertz;
                unit volt(unit_data(2, -3, 1, -1, 0, 0, 0));
                unit V = volt;
                unit newton(unit_data(1, -2, 1, 0, 0, 0, 0));
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

            } // namespace derived_units
            

            unit radians = unit(180 / math::constants::pi, special::one); 
            unit rad = radians;
            unit mps(m / s);
            unit mpss(m / s.pow(2)); 


            // unit hr(60.0, min, "hr");
            // unit day(24.0, hr, "day");
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

            using namespace physics::units; 

            // generate a conversion factor between two units in a constexpr function, the
            // units will only convert if they have the same base unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(UX start, UX2 result) { return quick_convert(1.0, start, result); }

            // generate a conversion factor between two units in a constexpr function, the units will only convert if they have the same base unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(double val, const UX& start, const UX2& result) {
                static_assert(std::is_same<UX, unit>::value || std::is_same<UX, unit>::value,
                    "convert argument types must be unit");
                static_assert(std::is_same<UX2, unit>::value || std::is_same<UX2, unit>::value,
                    "convert argument types must be unit");
                return (start.data() == result.data()) ?  val * start.multiplier() / result.multiplier() : 
                    constants::invalid_conversion;
            }

            // convert a value from one unit base to another
            template<typename UX, typename UX2>
            double convert(double val, const UX& start, const UX2& result) {
                static_assert(std::is_same<UX, unit>::value || std::is_same<UX, unit>::value,
                    "convert argument types must be unit");
                static_assert(std::is_same<UX2, unit>::value || std::is_same<UX2, unit>::value, 
                    "convert argument types must be unit");     
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
            double convert(const UX& start, const UX2& result) { return convert(1.0, start, result); }


        } // namespace tools


    } // namespace math


    namespace physics {

        // namespace defining some usefull measurement structures
        namespace measurements {

            using namespace units; 

            // measurement class with a value and an unit
            class measurement {

                protected:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{};

                    unit unit_;


                public:

                    // =============================================                                                                                         
                    // constructors & destructor 
                    // =============================================  

                    // default constructor
                    measurement() noexcept {};

                    // constructor from a value and a unit 
                    explicit constexpr measurement(const double& value, const unit& unit) noexcept : value_{value}, unit_{unit} {}

                    // constructor from a measurement
                    constexpr measurement(const measurement& other) : value_{other.value_}, unit_{other.unit_} {}
                        
                    // destructor
                    ~measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    constexpr measurement operator*(const measurement& other) const { return measurement(value_ * ((unit_ == other.unit_) ? other.value_ : other.value_as(unit_)), unit_ * other.unit_); }
                
                    constexpr measurement operator*(const double& val) const { return measurement(value_ * val, unit_); }
                    
                    constexpr measurement operator/(const measurement& other) const { return measurement(value_ / other.value_, unit_ / other.unit_); }

                    constexpr measurement operator/(const double& val) const { return measurement(value_ / val, unit_); } 

                    constexpr measurement& operator=(const measurement& other) noexcept {
                        value_ = (unit_ == other.unit_) ? other.value_ : other.value_as(unit_);
                        return *this;
                    }

                    constexpr measurement& operator=(const double& val) noexcept {
                        value_ = val;
                        return *this;
                    }     
    
                    constexpr bool operator==(const double& val) const { return (value_ == val) ? true : math::tools::compare_round_equals(value_, val); }

                    inline bool operator==(const measurement& other) const { return value_equality_check((unit_ == other.units()) ? other.value() : other.value_as(unit_)); }

                    constexpr bool operator!=(const double& val) const { return !operator==(val); }

                    inline bool operator!=(const measurement& other) const { return !value_equality_check((unit_ == other.units()) ? other.value() : other.value_as(unit_)); }

                    constexpr bool operator>(const double& val) const { return value_ > val; }

                    inline bool operator>(const measurement& other) const { return value_ > other.value_as(unit_); }

                    constexpr bool operator<(const double& val) const { return value_ < val; }

                    inline bool operator<(const measurement& other) const { return value_ < other.value_as(unit_); }
    
                    constexpr bool operator>=(const double& val) const { return (value_ >= val) ? true : operator==(val); }

                    inline bool operator>=(const measurement& other) const { return (value_ > other.value_as(unit_)) ? true : value_equality_check(other.value_as(unit_)); }

                    constexpr bool operator<=(const double& val) const { return value_ <= val ? true : operator==(val); }

                    inline bool operator<=(const measurement& other) const { return (value_ < other.value_as(unit_)) ? true : value_equality_check(other.value_as(unit_)); }       

                private:

                    // does a numerical equality check on the value accounting for tolerances
                    constexpr bool value_equality_check(const double& val) const {
                        return (value_ == val) ? true : math::tools::compare_round_equals(value_, val);
                    } 


                public:

                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a measurement 
                    constexpr measurement inv() const { return measurement(1 / value_, unit_.inv()); }

                    // take a measurement to a power
                    constexpr measurement pow(const int& power) const { return measurement(std::pow(value_, power), unit_.pow(power)); }

                    constexpr measurement square() const { return measurement(math::algebra::square(value_), unit_.square()); }

                    constexpr measurement cube() const { return measurement(math::algebra::cube(value_), unit_.cube()); }

                    // take a measurement to a root power
                    inline measurement root(const int& power) const { return measurement(math::algebra::root(value_, power), unit_.root(power)); }

                    inline measurement sqrt() const { return measurement(math::algebra::sqrt(value_), unit_.sqrt()); }

                    inline measurement cbrt() const { return measurement(math::algebra::cbrt(value_), unit_.cbrt()); }

                    
                    // =============================================                                                                                         
                    // set & get, convert & print methods
                    // =============================================  

                    // set the numerical value of the measurement
                    constexpr void value(const double& value) { value_ = value; }

                    // get the numerical value of the measurement
                    constexpr double value() const { return value_; }                        
                    
                    // get the numerical value as a particular unit type
                    inline double value_as(const unit& desired_units) const { return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units); }

                    // get the unit component of a measurement
                    constexpr unit units() const { return unit_; }

                    // convert the measurement to a single unit
                    constexpr unit as_unit() const { return unit(value_, unit_); }

                    // convert a unit to have a new base
                    inline measurement convert_to(const unit& newUnits) const { return measurement(math::tools::convert(value_, unit_, newUnits), newUnits); }

                    // convert a unit into its base units
                    constexpr measurement convert_to_base() const { return measurement(value_ * unit_.multiplier(), unit_.as_unit()); }

                    // print the measurement
                    void print() const { 
                        // if (unit_.multiplier()!= 1) std::cout << std::setprecision(unit_.multiplier()) << value_ << " "; 
                        std::cout << value_ << " ";
                        unit_.print(); 
                        std::cout << "\n"; 
                    }


            }; // class measurement


            // measurement class with a value and a fixed unit
            class fixed_measurement {

                protected:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{};  

                    const unit unit_;  


                public:

                    // =============================================                                                                                         
                    // constructors & destructor
                    // =============================================  

                    // default constructor
                    fixed_measurement() noexcept {};

                    // constructor from a value and a unit 
                    explicit constexpr fixed_measurement(const double& val, const unit& unit) noexcept : value_(val), unit_(unit) {}

                    // constructor from a measurement
                    constexpr fixed_measurement(const measurement& other) : value_(other.value()), unit_(other.units()) {}

                    // constructor from a fixed measurement
                    constexpr fixed_measurement(const fixed_measurement& other) : value_(other.value_), unit_(other.unit_) {}

                    // destructor
                    ~fixed_measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    constexpr fixed_measurement operator+(const double& val) const {
                        return fixed_measurement(value_ + val, unit_);
                    }
                    
                    constexpr fixed_measurement operator-(const double& val) const {
                        return fixed_measurement(value_ - val, unit_);
                    }

                    constexpr fixed_measurement& operator+=(const double& val) {
                        value_ += val;
                        return *this;
                    }

                    constexpr fixed_measurement& operator-=(const double& val) {
                        value_ -= val;
                        return *this;
                    }
                    
                    constexpr fixed_measurement operator*(const double& val) const {
                        return fixed_measurement(value_ * val, unit_);
                    }

                    constexpr fixed_measurement operator/(const double& val) const {
                        return fixed_measurement(value_ / val, unit_);
                    }
    
                    constexpr fixed_measurement& operator*=(const double& val) {
                        value_ *= val;
                        return *this;
                    }

                    constexpr fixed_measurement& operator/=(const double& val) {
                        value_ /= val;
                        return *this;
                    }

                    inline fixed_measurement operator+(const fixed_measurement& other) const {
                        return fixed_measurement(value_ + other.value_as(unit_), unit_);
                    }

                    inline fixed_measurement operator+(const measurement& other) const {
                        return fixed_measurement(value_ + other.value_as(unit_), unit_);
                    }

                    inline fixed_measurement operator-(const fixed_measurement& other) const {
                        return fixed_measurement(value_ - other.value_as(unit_), unit_);
                    }

                    inline fixed_measurement operator-(const measurement& other) const {
                        return fixed_measurement(value_ - other.value_as(unit_), unit_);
                    }   

                    constexpr fixed_measurement operator*(const fixed_measurement& other) const {
                        return fixed_measurement(value_ * other.value_, unit_ * other.unit_);
                    }

                    constexpr fixed_measurement operator/(const fixed_measurement& other) const {
                        return fixed_measurement(value_ / other.value_, unit_ / other.unit_);
                    }

                    constexpr fixed_measurement& operator=(const fixed_measurement& other) noexcept {
                        value_ = (unit_ == other.units()) ? other.value() : other.value_as(unit_);
                        return *this;
                    }

                    constexpr fixed_measurement& operator=(const measurement& other) noexcept {
                        value_ = (unit_ == other.units()) ? other.value() : other.value_as(unit_);
                        return *this;
                    }

                    constexpr fixed_measurement& operator=(const double& val) noexcept {
                        value_ = val;
                        return *this;
                    }  

                    inline bool operator==(const fixed_measurement& val) const {
                        return operator==((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator==(const measurement& val) const {
                        return operator==((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }
                    
                    constexpr bool operator==(const double& val) const {
                        return (value_ == val) ? true : math::tools::compare_round_equals(value_, val);
                    }
                    
                    inline bool operator!=(const fixed_measurement& val) const {
                        return operator!=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator!=(const measurement& val) const {
                        return operator!=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    constexpr bool operator!=(const double& val) const { return !operator==(val); }

                    inline bool operator>(const fixed_measurement& val) const {
                        return operator>((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator>(const measurement& val) const {
                        return operator>((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    constexpr bool operator>(const double& val) const { return value_ > val; }

                    inline bool operator<(const fixed_measurement& val) const {
                        return operator<((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator<(const measurement& val) const {
                        return operator<((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    constexpr bool operator<(const double& val) const { return value_ > val; }

                    inline bool operator>=(const fixed_measurement& val) const {
                        return operator>=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator>=(const measurement& val) const {
                        return operator>=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    constexpr bool operator>=(const double& val) const { return (value_ >= val) ? true : operator==(val); }

                    inline bool operator<=(const fixed_measurement& val) const {
                        return operator<=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator<=(const measurement& val) const {
                        return operator<=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    constexpr bool operator<=(const double& val) const { return (value_ <= val) ? true : operator==(val); }


                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a fixed_measurement 
                    constexpr fixed_measurement inv() const { return fixed_measurement(1 / value_, unit_.inv()); }

                    // take a fixed_measurement to a power
                    constexpr fixed_measurement pow(const int& power) const { return fixed_measurement(std::pow(value_, power), unit_.pow(power)); }

                    constexpr fixed_measurement square() const { return fixed_measurement(math::algebra::square(value_), unit_.square()); }

                    constexpr fixed_measurement cube() const { return fixed_measurement(math::algebra::cube(value_), unit_.cube()); }

                    // take a fixed_measurement to a root power
                    inline fixed_measurement root(const int& power) const { return fixed_measurement(math::algebra::root(value_, power), unit_.root(power)); }

                    inline fixed_measurement sqrt() const { return fixed_measurement(math::algebra::sqrt(value_), unit_.sqrt()); }

                    inline fixed_measurement cbrt() const { return fixed_measurement(math::algebra::cbrt(value_), unit_.cbrt()); }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // direct conversion operator
                    constexpr operator measurement() { return measurement(value_, unit_); }

                    // convert a unit to have a new base
                    inline measurement convert_to(const unit& newUnits) const {
                        return measurement(math::tools::convert(value_, unit_, newUnits), newUnits);
                    }

                    // convert a unit into its base units
                    inline measurement convert_to_base() const {
                        return measurement(value_ * unit_.multiplier(), unit(unit_.data()));
                    }


                    // =============================================                                                                                         
                    // set & get & print methods
                    // =============================================  

                    // set the numerical component of the measurement
                    constexpr void value(const double& value) { value_ = value; }

                    // get the numerical component of the measurement
                    constexpr double value() const { return value_; }                        
                    
                    // get the numerical value as a particular unit type
                    double value_as(const unit& desired_units) const {
                        return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units);
                    }

                    // get the unit component of a measurement
                    constexpr unit units() const { return unit_; }

                    // convert the measurement to a single unit
                    constexpr unit as_unit() const { return {value_, unit_ }; }

                    // print the measurement
                    void print() const { 
                        // if (unit_.multiplier() != 1) std::cout << std::setprecision(unit_.multiplier()) << value_ << " "; 
                        std::cout << value_ << " ";
                        unit_.print(); 
                        std::cout << "\n"; 
                    }


            }; // class fixed_measurement


            // measurement class with an uncertain value and an unit
            class uncertain_measurement {

                private:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{}, uncertainty_{};

                    unit unit_;       
                

                public: 

                    // =============================================                                                                                         
                    // constructors & destructor 
                    // =============================================  

                    // default constructor
                    uncertain_measurement() noexcept {};

                    // constructor from a value, uncertainty, and unit
                    explicit constexpr uncertain_measurement(const double& val, const double& uncertainty_val, const unit& unit) noexcept :
                        value_(val), uncertainty_(uncertainty_val), unit_(unit) {}

                    // constructpr from a value and an unit, assuming the uncertainty is 0
                    explicit constexpr uncertain_measurement(const double& val, const unit& unit) noexcept :
                        value_(val), unit_(unit) {}

                    // constructor from a measurement and uncertainty value
                    explicit constexpr uncertain_measurement(const measurement& other, const double& uncertainty_val) noexcept : 
                        value_(other.value()), uncertainty_(uncertainty_val), unit_(other.units()) {}

                    // constructor from a fixed_measurement and uncertainty value
                    explicit constexpr uncertain_measurement(const fixed_measurement& other, const double& uncertainty_val) noexcept : 
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
                    uncertain_measurement operator*(const uncertain_measurement& other) const {
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

                    constexpr uncertain_measurement operator*(const unit& other) const {
                        return uncertain_measurement(value_, uncertainty_, unit_ * other);
                    }

                    constexpr uncertain_measurement operator*(const double& val) const {
                        return uncertain_measurement(value_ * val, uncertainty_ * val, unit_);
                    }

                    // compute a unit division and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator/(const uncertain_measurement& other) const {
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

                    constexpr uncertain_measurement operator/(const unit& other) const {
                        return uncertain_measurement(value_, uncertainty_, unit_ / other);
                    }

                    constexpr uncertain_measurement operator/(const double& val) const {
                        return uncertain_measurement(value_ / val, uncertainty_ / val, unit_);
                    }

                    // compute a unit addition and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator+(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                        double ntol = math::algebra::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                        return uncertain_measurement(value_ + cval * other.value_, ntol, unit_);
                    }

                    // compute a unit addition and calculate the new uncertainties using the simple uncertainty summation method
                    uncertain_measurement simple_add(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return uncertain_measurement(value_ + cval * other.value_, ntol, unit_);
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator-(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                        double ntol = math::algebra::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                        return uncertain_measurement(value_ - cval * other.value_, ntol, unit_);
                    }

                    // compute a unit subtraction and calculate the new uncertainties using the simple uncertainty summation method
                    uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                        auto cval = math::tools::convert(other.unit_, unit_);
                        double ntol = uncertainty_ + other.uncertainty_ * cval;
                        return uncertain_measurement(value_ - cval * other.value_, ntol, unit_);
                    }

                    inline uncertain_measurement operator+(const measurement& other) const {
                        return uncertain_measurement(value_ + other.value_as(unit_), uncertainty_, unit_);
                    }

                    inline uncertain_measurement operator-(const measurement& other) const {
                        return uncertain_measurement(value_ - other.value_as(unit_), uncertainty_, unit_);
                    }

                    // comparison operators 
                    bool operator==(const measurement& other) const {
                        if (uncertainty_ == 0.0F) { return (value_ == other.value_as(unit_)) ? true : math::tools::compare_round_equals(value_, other.value_as(unit_)); }
                        return (other.value_as(unit_) >= (value_ - uncertainty_) && other.value_as(unit_) <= (value_ + uncertainty_));
                    }

                    inline bool operator>(const measurement& other) const { return value_ > other.value_as(unit_); }

                    inline bool operator<(const measurement& other) const { return value_ < other.value_as(unit_); }

                    inline bool operator>=(const measurement& other) const {
                        return (value() >= other.value_as(unit_)) ? true : operator==(measurement(other.value_as(unit_), unit_));
                    }

                    inline bool operator<=(const measurement& other) const {
                        return (value() <= other.value_as(unit_)) ? true : operator==(measurement(other.value_as(unit_), unit_));
                    }

                    inline bool operator!=(const measurement& other) const { return !operator==(other); }

                    inline bool operator==(const uncertain_measurement& other) const { return (simple_subtract(other) == measurement(0.0, unit_)); }

                    inline bool operator>(const uncertain_measurement& other) const { return value_ > other.value_as(unit_); }

                    inline bool operator<(const uncertain_measurement& other) const { return value_ < other.value_as(unit_); }

                    inline bool operator>=(const uncertain_measurement& other) const {
                        return (simple_subtract(other).value_ >= 0.0F) ? true : (simple_subtract(other) == measurement(0.0, unit_));
                    }

                    inline bool operator<=(const uncertain_measurement& other) const {
                        return (simple_subtract(other).value_ <= 0.0F) ? true : (simple_subtract(other) == measurement(0.0, unit_));
                    }

                    inline bool operator!=(const uncertain_measurement& other) const { return !operator==(other); }


                    // =============================================                                                                                         
                    // get and set methods
                    // =============================================  
                    
                    // set the uncertainty
                    constexpr inline void uncertainty(const double& uncert) { uncertainty_ = uncert; }

                    // set the uncertainty
                    constexpr uncertain_measurement& uncertainty(double newUncertainty) {
                        uncertainty_ = newUncertainty;
                        return *this;
                    }

                    // set the uncertainty
                    uncertain_measurement& uncertainty(const measurement& newUncertainty) {
                        uncertainty_ = newUncertainty.value_as(unit_);
                        return *this;
                    }

                    // get the uncertainty as a separate measurement
                    constexpr measurement uncertainty_measurement() const { return uncertain_measurement(uncertainty(), unit_); }

                    // // cast operator to a measurement
                    constexpr operator measurement() const { return measurement(value(), unit_); }

                    // get the numerical value 
                    constexpr double value() const { return value_; }

                    // get the numerical value of the uncertainty
                    constexpr double uncertainty() const { return uncertainty_; }

                    // get the underlying units value
                    constexpr unit units() const { return unit_; }

                    // get the numerical value as a particular unit type
                    inline double value_as(const unit& desired_units) const {
                        return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units);
                    }

                    // get the numerical value of the uncertainty as a particular unit
                    inline double uncertainty_as(const unit& desired_units) const {
                        return (unit_ == desired_units) ? uncertainty_ : math::tools::convert(uncertainty_, unit_, desired_units);
                    }

                    // convert a unit to have a new base
                    uncertain_measurement convert_to(const unit& newUnits) const {
                        auto cval = math::tools::convert(unit_, newUnits);
                        return uncertain_measurement(cval * value_, uncertainty_ * cval, newUnits);
                    }

                    // print the uncertain measurement
                    void print() const { 
                        // if (unit_.multiplier() != 1) std::cout << std::setprecision(unit_.multiplier()) << value() << " ± " << uncertainty() << " "; 
                        std::cout << value() << " ± " << uncertainty() << " "; 
                        unit_.print(); 
                        std::cout << "\n";
                    }


            }; // class uncertain_measurement


        } // namespace measurements


    } // namespace physics


    namespace math {

        // namespace defining some usefull operation
        namespace algebra {

            using namespace physics::units; 
            using namespace physics::measurements;

            // generate a unit which is an integer power of another
            constexpr unit pow(const unit& u, const int& power) { return u.pow(power); }

            // generate the square of an unit
            constexpr unit square(const unit& u) { return u.pow(2); }

            // generate the cube of an unit
            constexpr unit cube(const unit& u) { return u.pow(3); }

            // generate the root of an unit
            unit root(const unit& un, const int& power) {
                if (power == 0) { return special::one; }
                if (un.multiplier() < 0.0 && power % 2 == 0) { return special::invalid; }
                return unit(un.data().root(power), root(un.multiplier(), power));
            }

            // generate the square root of an unit
            inline unit sqrt(const unit& u) { return sqrt(u); }

            // generate the cubic root of an unit
            inline unit cbrt(const unit& u) { return cbrt(u); }

            constexpr measurement pow(const measurement& meas, const int& power) {
                return measurement{ std::pow(meas.value(), power), meas.units().pow(power) };
            }

            constexpr fixed_measurement pow(const fixed_measurement& meas, const int& power) {
                return fixed_measurement{std::pow(meas.value(), power), meas.units().pow(power)};
            }

            uncertain_measurement pow(const uncertain_measurement& meas, const int& power) {
                auto new_value = std::pow(meas.value(), power);
                auto new_tol = ((power >= 0) ? power : -power) * new_value * meas.uncertainty() / meas.value();
                return uncertain_measurement(new_value, new_tol, meas.units().pow(power));                    
            }

            constexpr measurement square(const measurement& meas) { return pow(meas, 2); }
            
            constexpr fixed_measurement square(const fixed_measurement& meas) { return pow(meas, 2); }

            inline uncertain_measurement square(const uncertain_measurement& meas) { return pow(meas, 2); }

            constexpr measurement cube(const measurement& meas) { return pow(meas, 3); }
            
            constexpr fixed_measurement cube(const fixed_measurement& meas) { return pow(meas, 3); }
            
            inline uncertain_measurement cube(const uncertain_measurement& meas) { return pow(meas, 3); }


            inline measurement root(const measurement& meas, const int& power) {
                return measurement(root(meas.value(), power), root(meas.units(), power));
            }

            inline fixed_measurement root(const fixed_measurement& meas, const int& power) {
                return fixed_measurement(root(meas.value(), power), root(meas.units(), power));
            }

            uncertain_measurement root(const uncertain_measurement& um, const int& power) {
                auto new_value = root(um.value(), power);
                auto new_tol = new_value * um.uncertainty() / (static_cast<double>((power >= 0) ? power : -power) * um.value());
                return uncertain_measurement(new_value, new_tol, root(um.units(), power));
            }

            inline measurement sqrt(const measurement& meas) { return root(meas, 2); }
            
            inline fixed_measurement sqrt(const fixed_measurement& meas) { return root(meas, 2); }

            inline uncertain_measurement sqrt(const uncertain_measurement& meas) { return root(meas, 2); }            

            inline measurement cbrt(const measurement& meas) { return root(meas, 3); }

            inline fixed_measurement cbrt(const fixed_measurement& meas) { return root(meas, 3); }

            inline uncertain_measurement cbrt(const uncertain_measurement& meas) { return root(meas, 3); }


            constexpr bool operator==(const double& val, const fixed_measurement& v2) { return v2 == val; }

            constexpr bool operator!=(const double& val, const fixed_measurement& v2) { return v2 != val; }

            constexpr bool operator>(const double& val, const fixed_measurement& v2) { return val > v2.value(); }

            constexpr bool operator<(const double& val, const fixed_measurement& v2) { return val < v2.value(); }

            constexpr bool operator>=(const double& val, const fixed_measurement& v2) { return (val >= v2.value()) ? true : (v2 == val); }

            constexpr bool operator<=(const double& val, const fixed_measurement& v2) { return (val <= v2.value()) ? true : (v2 == val); }

            inline bool operator==(const measurement& other, const uncertain_measurement& v2) { return v2 == other; }

            inline bool operator!=(const measurement& other, const uncertain_measurement& v2) { return v2 != other; }

            constexpr bool operator>(const measurement& other, const uncertain_measurement& v2) { return other.value() > v2.value(); }

            constexpr bool operator<(const measurement& other, const uncertain_measurement& v2) { return other.value() < v2.value(); }

            constexpr bool operator>=(const measurement& other, const uncertain_measurement& v2) { return (other > v2) ? true : (v2 == other); }

            constexpr bool operator<=(const measurement& other, const uncertain_measurement& v2) { return (other < v2) ? true : (v2 == other); }

    
        } // namespace algebra


    } // namespace math


            // constexpr physics::measurements::measurement operator*(const double& val, const physics::measurements::measurement& meas) { return meas * val; }

            // constexpr physics::measurements::measurement operator*(const double& val, const physics::units::unit& unit_base) { return physics::measurements::measurement(val, unit_base); }

            // constexpr physics::measurements::measurement operator*(const physics::units::unit& unit_base, const double& val) { return physics::measurements::measurement(val, unit_base); }

            // constexpr physics::measurements::fixed_measurement operator*(const double& v1, const physics::measurements::fixed_measurement& v2) { return physics::measurements::fixed_measurement(v1 * v2.value(), v2.units()); }

            // constexpr physics::measurements::fixed_measurement operator*(const physics::measurements::fixed_measurement& v2, const double& v1) { return physics::measurements::fixed_measurement(v1 * v2.value(), v2.units()); }

            // inline physics::measurements::uncertain_measurement operator*(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) { return v2.operator*(v1); }

            // inline physics::measurements::uncertain_measurement operator*(const double& v1, const physics::measurements::uncertain_measurement& v2) { return v2.operator*(v1); }

            // inline physics::measurements::uncertain_measurement operator*(const physics::measurements::uncertain_measurement& v2, const double& v1) { return v2.operator*(v1); }


            // constexpr physics::measurements::measurement operator/(const double& val, const physics::measurements::measurement& meas) { return physics::measurements::measurement(val / meas.value(), meas.units().inv()); }

            // constexpr physics::measurements::measurement operator/(const physics::measurements::measurement& meas, const double& val) { return physics::measurements::measurement(val / meas.value(), meas.units().inv()); }
            
            // constexpr physics::measurements::measurement operator/(const double& val, const physics::units::unit& unit_base) { return physics::measurements::measurement(val, unit_base.inv()); }
            
            // constexpr physics::measurements::measurement operator/(const physics::units::unit& unit_base, const double& val) { return physics::measurements::measurement(1.0 / val, unit_base); }

            // constexpr physics::measurements::fixed_measurement operator/(const double& v1, const physics::measurements::fixed_measurement& v2) { return physics::measurements::fixed_measurement(v1 / v2.value(), v2.units().inv()); }             
            
            // constexpr physics::measurements::fixed_measurement operator/(const physics::measurements::fixed_measurement& v2, const double& v1) { return physics::measurements::fixed_measurement(v1 / v2.value(), v2.units().inv()); }

            // physics::measurements::uncertain_measurement operator/(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
            //     double ntol = v2.uncertainty() / v2.value();
            //     double nval = v1.value() / v2.value();
            //     return physics::measurements::uncertain_measurement(nval, nval * ntol, v1.units() / v2.units());
            // }

            // physics::measurements::uncertain_measurement operator/(const double& v1, const physics::measurements::uncertain_measurement& v2) {
            //     double ntol = v2.uncertainty() / v2.value();
            //     double nval = v1 / v2.value();
            //     return physics::measurements::uncertain_measurement(nval, nval * ntol, v2.units().inv());
            // }

            // physics::measurements::uncertain_measurement operator/(const int& v1, const physics::measurements::uncertain_measurement& v2) {
            //     double ntol = v2.uncertainty() / v2.value();
            //     double nval = static_cast<double>(v1) / v2.value();
            //     return physics::measurements::uncertain_measurement(nval, nval * ntol, v2.units().inv());
            // }


            // constexpr physics::measurements::fixed_measurement operator+(const double& v1, const physics::measurements::fixed_measurement& v2) { return physics::measurements::fixed_measurement(v1 + v2.value(), v2.units()); }

            // constexpr physics::measurements::fixed_measurement operator-(const double& v1, const physics::measurements::fixed_measurement& v2) { return physics::measurements::fixed_measurement(v1 - v2.value(), v2.units()); }
            
            // physics::measurements::uncertain_measurement operator+(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
            //     double cval = math::tools::convert(v2.units(), v1.units());
            //     double ntol = v2.uncertainty() * cval;
            //     return physics::measurements::uncertain_measurement(v1.value() + cval * v2.value(), ntol, v1.units());
            // }

            // physics::measurements::uncertain_measurement operator-(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
            //     double cval = math::tools::convert(v2.units(), v1.units());
            //     double ntol = v2.uncertainty() * cval;
            //     return physics::measurements::uncertain_measurement(v1.value() - cval * v2.value(), ntol, v1.units());
            // }


    // operations amongs std::vector<T>
    template <typename T> 
    std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] + v2[i];
        return vec;
    }

    template <typename T> 
    std::vector<T> operator+(const std::vector<T>& v1, const double& val) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] + val;
        return vec;
    }

    template <typename T> 
    std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] - v2[i];
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator-(const std::vector<T>& v1, const double& val) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] - val;
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator+=(const std::vector<T>& v1, const std::vector<T>& v2) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] + v2[i];
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator+=(const std::vector<T>& v1, const double& val) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] + val;
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator-=(const std::vector<T>& v1, const std::vector<T>& v2) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] - v2[i];
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator-=(const std::vector<T>& v1, const double& val) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] - val;
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T>& v2) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] * v2[i];
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator*(const std::vector<T>& v1, const double& val) {
        std::cout << "operator* called\n"; 
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = (T) (v1[i] * val);
        return vec;
    }
    
    template <typename T>            
    std::vector<T> operator*(const double& val, const std::vector<T>& v1) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] * val;
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator/(const std::vector<T>& v1, const std::vector<T>& v2) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] / v2[i];
        return vec;
    }
    
    template <typename T> 
    std::vector<T> operator/(const std::vector<T>& v1, const double& val) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] / val;
        return vec;
    }
    
    template <typename T>            
    std::vector<T> operator/(const double& val, const std::vector<T>& v1) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = v1[i] / val;
        return vec;
    }
    
    template <typename T>            
    std::vector<T> pow(const std::vector<T>& v1, const double& power) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = math::algebra::pow(v1[i], power);
        return vec;
    }
    
    template <typename T>            
    std::vector<T> square(const std::vector<T>& v1) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = math::algebra::square(v1[i]);
        return vec;
    }
    
    template <typename T>            
    std::vector<T> cube(const std::vector<T>& v1) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = math::algebra::cube(v1[i]);
        return vec;
    }
    
    template <typename T>            
    std::vector<T> root(const std::vector<T>& v1, const double& power) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = math::algebra::root(v1[i], power);
        return vec;
    }
    
    template <typename T>            
    std::vector<T> sqrt(const std::vector<T>& v1) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = math::algebra::sqrt(v1[i]);
        return vec;
    }
    
    template <typename T>            
    std::vector<T> cbrt(const std::vector<T>& v1) {
        std::vector<T> vec = std::vector<T>(v1.size());
        for (size_t i{}; i < v1.size(); i++) vec[i] = math::algebra::cbrt(v1[i]);
        return vec;
    }


    // operations amongs std::pair<T1, T2>
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator+(const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
        return std::make_pair((T1) v1.first + v2.first, (T2) v1.second + v2.second);
    }

    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator+(const std::pair<T1, T2>& v1, const double& val) {
        return std::make_pair(v1.first + val, v1.second + val);
    }

    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator-(const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
        return std::make_pair(v1.first - v2.first, v1.second - v2.second);
    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator-(const std::pair<T1, T2>& v1, const double& val) {
        return std::make_pair(v1.first - val, v1.second - val);
    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator+=(const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
        return std::make_pair(v1.first + v2.first, v1.second + v2.second);

    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator+=(const std::pair<T1, T2>& v1, const double& val) {
        return std::make_pair(v1.first + val, v1.second + val);
    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator-=(const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
        return std::make_pair(v1.first + v2.first, v1.second + v2.second);
    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator-=(const std::pair<T1, T2>& v1, const double& val) {
        return std::make_pair(v1.first + val, v1.second + val);
    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator*(const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
        return std::make_pair(static_cast<T1>(v1.first * v2.first)), static_cast<T2>(v1.second * v2.second);
    }

    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator*(const std::pair<T1, T2>& v1, const double& val) {
        return std::make_pair(static_cast<T1>(v1.first) * val, static_cast<T2>(v1.second) * val);
    }

    template <typename T1, typename T2>            
    inline std::pair<T1, T2> operator*(const double& val, const std::pair<T1, T2>& v1) {
        return std::make_pair(static_cast<T1>(v1.first) * val, static_cast<T2>(v1.second) * val);
    }
    
    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator/(const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2) {
        return std::make_pair(v1.first / v2.first, v1.second / v2.second);
    }

    template <typename T1, typename T2> 
    inline std::pair<T1, T2> operator/(const std::pair<T1, T2>& v1, const double& val) {
        return std::make_pair(static_cast<T1>(v1.first) / val, static_cast<T2>(v1.second) / val);
    }

    template <typename T1, typename T2>            
    inline std::pair<T1, T2> operator/(const double& val, const std::pair<T1, T2>& v1) {
        return std::make_pair(v1.first / val, v1.second / val);
    }
    
    template <typename T1, typename T2>            
    inline std::pair<T1, T2> pow(const std::pair<T1, T2>& v1, const double& power) {
        return std::make_pair(math::algebra::pow(v1.first, power), math::algebra::pow(v1.secon, power));
    }
    
    template <typename T1, typename T2>            
    inline std::pair<T1, T2> square(const std::pair<T1, T2>& v1) {
        return std::make_pair(math::algebra::square(v1.first), math::algebra::square(v1.second));
    }
    
    template <typename T1, typename T2>            
    inline std::pair<T1, T2> cube(const std::pair<T1, T2>& v1) {
        return std::make_pair(math::algebra::cube(v1.first), math::algebra::cube(v1.second));
    }
    
    template <typename T1, typename T2>            
    inline std::pair<T1, T2> root(const std::pair<T1, T2>& v1, const double& power) {
        return std::make_pair(math::algebra::root(v1.first, power), math::algebra::root(v1.secon, power));
    }
    
    template <typename T1, typename T2>            
    inline std::pair<T1, T2> sqrt(const std::pair<T1, T2>& v1) {
        return std::make_pair(math::algebra::sqrt(v1.first), math::algebra::sqrt(v1.second));
    }
    
    template <typename T1, typename T2>            
    inline std::pair<T1, T2> cbrt(const std::pair<T1, T2>& v1) {
        return std::make_pair(math::algebra::cbrt(v1.first), math::algebra::cbrt(v1.second));
    }
    

    namespace physics {


        // class expressing the coordinate of an object in a 3D system
        class position {

            private: 
                
                // =============================================
                // class members
                // =============================================
                
                std::vector<measurements::fixed_measurement> pos_ = std::vector<measurements::fixed_measurement>(3);


            public: 

                // =============================================
                // constructors & destructor
                // =============================================
                
                position() noexcept {}

                explicit position(const measurements::fixed_measurement& x, const measurements::fixed_measurement& y, const measurements::fixed_measurement& z) noexcept : 
                    pos_{x, y, z} {}

                explicit position(const double& x, const double& y, const double& z) noexcept : 
                    pos_{measurements::fixed_measurement(x, units::m), measurements::fixed_measurement(y, units::m), measurements::fixed_measurement(z, units::m)} {}

                position(const std::vector<measurements::fixed_measurement>& pos) noexcept : pos_{pos} {}

                position(const position& pos) noexcept : position(pos.get_position()) {}

                ~position() = default; 


                // =============================================
                // operators
                // =============================================

                inline position operator+(const position& other) { return position(pos_ + other.pos_); } 

                inline position operator+(const double& val) { return position(pos_ + val); }

                inline position operator-(const position& other) { return position(pos_ - other.pos_); } 

                inline position operator-(const double& val) { return position(pos_ - val); }

                inline position operator*(const position& other) { return position(pos_ * other.pos_); } 

                inline position operator*(const double& val) { return position(pos_ * val); }

                inline position operator/(const position& other) { return position(pos_ / other.pos_); } 

                inline position operator/(const double& val) { return position(pos_ / val); }


                // =============================================
                // operations
                // =============================================


                // =============================================
                // set, get and print methods
                // =============================================

                inline void set_position(const measurements::fixed_measurement& x, const measurements::fixed_measurement& y, const measurements::fixed_measurement& z) { 
                    pos_ = std::vector<measurements::fixed_measurement>{x, y, z}; 
                }

                inline void set_position(const double& x, const double& y, const double& z) { 
                    pos_ = std::vector<measurements::fixed_measurement>{measurements::fixed_measurement(x, units::m), measurements::fixed_measurement(y, units::m), measurements::fixed_measurement(z, units::m)}; 
                }

                inline void set_position(const std::vector<measurements::fixed_measurement>& pos) { pos_ = pos; }
                
                inline void set_position(const position& pos) { pos_ = pos.get_position(); }

                inline void x(const measurements::fixed_measurement& x) { pos_[0] = x; }
                
                inline void y(const measurements::fixed_measurement& y) { pos_[1] = y; }
                
                inline void z(const measurements::fixed_measurement& z) { pos_[2] = z; }

                inline void x(const double& x) { pos_[0] = measurements::fixed_measurement(x, pos_[0].units()); }

                inline void y(const double& y) { pos_[1] = measurements::fixed_measurement(y, pos_[1].units()); }

                inline void z(const double& z) { pos_[2] = measurements::fixed_measurement(z, pos_[2].units()); }

                inline std::vector<measurements::fixed_measurement> get_position() const { return pos_; }

                inline measurements::fixed_measurement get_x() const { return pos_[0]; }
                
                inline measurements::fixed_measurement get_y() const { return pos_[1]; }
                
                inline measurements::fixed_measurement get_z() const { return pos_[2]; }             

                inline measurements::fixed_measurement magnitude() const { 
                    return math::algebra::sqrt(math::algebra::square(pos_[0]) + math::algebra::square(pos_[1]) + math::algebra::square(pos_[2]));
                }

                inline measurements::fixed_measurement distance(const position& other) const {        
                    return math::algebra::sqrt(math::algebra::square(other.pos_[0] - pos_[0]) + math::algebra::square(other.pos_[1] - pos_[1]) + math::algebra::square(other.pos_[2] - pos_[2]));
                }       

                inline measurements::fixed_measurement rho() const { 
                    return math::algebra::sqrt(math::algebra::square(pos_[0]) + math::algebra::square(pos_[1])); 
                }

                inline measurements::fixed_measurement phi() const { 
                    return measurements::fixed_measurement(atan2(pos_[1].value(), pos_[0].value()), units::rad); 
                }     

                inline measurements::fixed_measurement phi(const position& other) const { 
                    return measurements::fixed_measurement(atan2(other.pos_[1].value() - pos_[1].value(), other.pos_[0].value() - pos_[0].value()), units::rad); 
                }

                inline measurements::fixed_measurement theta() const { 
                    if (pos_[2] == 0) return measurements::fixed_measurement(0, units::rad);
                    return measurements::fixed_measurement(acos(pos_[2].value() / magnitude().value()), units::rad); 
                }

                inline measurements::fixed_measurement theta(const position& other) const {
                    if (other.pos_[2] == pos_[2]) return measurements::fixed_measurement(0, units::rad);
                    return measurements::fixed_measurement(acos((other.pos_[2].value() - pos_[2].value()) / distance(other).value()), units::rad); 
                }

                inline std::vector<double> direction() const {
                    return { cos(phi().value()), sin(phi().value()), pos_[2].value() / magnitude().value() };
                } 

                inline std::vector<double> direction(const position& other) const {
                    return { cos(phi(other).value()), sin(phi(other).value()), (other.pos_[2] - pos_[2]).value() / distance(other).value() };
                } 

                void print_position() const {
                    std::cout << "position = {\n";
                    for (auto i : pos_) i.print(); 
                    std::cout << "}\n"; 
                }


        }; // class position

        
        // class expressing the velocity of an object in a 3D system
        class velocity {

            private: 
                
                // =============================================
                // class members
                // =============================================
                
                std::vector<measurements::fixed_measurement> vel_ = std::vector<measurements::fixed_measurement>(3);
            

            public:  

                // =============================================
                // constructors
                // =============================================

                velocity() noexcept {}

                explicit velocity(const measurements::fixed_measurement& x, const measurements::fixed_measurement& y, const measurements::fixed_measurement& z) noexcept : 
                    vel_{x, y, z} {}

                explicit velocity(const double& x, const double& y, const double& z) noexcept : 
                    vel_{measurements::fixed_measurement(x, units::m), measurements::fixed_measurement(y, units::m), measurements::fixed_measurement(z, units::m)} {}

                velocity(const std::vector<measurements::fixed_measurement>& pos) : vel_{pos} {}

                velocity(const velocity& pos) : velocity(pos.get_velocity()) {}

                ~velocity() = default; 


                // =============================================
                // operators
                // =============================================

                inline velocity operator+(const velocity& other) { return velocity(vel_ + other.vel_); } 

                inline velocity operator+(const double& val) { return velocity(vel_ + val); }

                inline velocity operator-(const velocity& other) { return velocity(vel_ - other.vel_); } 

                inline velocity operator-(const double& val) { return velocity(vel_ - val); }

                inline velocity operator*(const velocity& other) { return velocity(vel_ * other.vel_); } 

                inline velocity operator*(const double& val) { return velocity(vel_ * val); }

                inline velocity operator/(const velocity& other) { return velocity(vel_ / other.vel_); } 

                inline velocity operator/(const double& val) { return velocity(vel_ / val); }


                // =============================================
                // set & get & print methods
                // =============================================

                inline void set_velocity(const measurements::fixed_measurement& x, const measurements::fixed_measurement& y, const measurements::fixed_measurement& z) { 
                    vel_ = std::vector<measurements::fixed_measurement>{x, y, z}; 
                }

                inline void set_velocity(const double& x, const double& y, const double& z) { 
                    vel_ = std::vector<measurements::fixed_measurement>{measurements::fixed_measurement(x, units::m), measurements::fixed_measurement(y, units::m), measurements::fixed_measurement(z, units::m)}; 
                }

                inline void set_velocity(const std::vector<measurements::fixed_measurement>& pos) { vel_ = pos; }
                
                inline void set_velocity(const velocity& pos) { vel_ = pos.get_velocity(); }

                inline void x(const measurements::fixed_measurement& x) { vel_[0] = x; }
                
                inline void y(const measurements::fixed_measurement& y) { vel_[1] = y; }
                
                inline void z(const measurements::fixed_measurement& z) { vel_[2] = z; }

                inline void x(const double& x) { vel_[0] = measurements::fixed_measurement(x, vel_[0].units()); }

                inline void y(const double& y) { vel_[1] = measurements::fixed_measurement(y, vel_[1].units()); }

                inline void z(const double& z) { vel_[2] = measurements::fixed_measurement(z, vel_[2].units()); }

                inline std::vector<measurements::fixed_measurement> get_velocity() const { return vel_; }

                inline measurements::fixed_measurement get_vx() const { return vel_[0]; }
                
                inline measurements::fixed_measurement get_vy() const { return vel_[1]; }
                
                inline measurements::fixed_measurement get_vz() const { return vel_[2]; }             

                inline measurements::fixed_measurement magnitude() const { 
                    return math::algebra::sqrt(math::algebra::square(vel_[0]) +
                                                math::algebra::square(vel_[1]) +
                                                math::algebra::square(vel_[2]));
                }
 
                inline measurements::fixed_measurement phi() const { 
                    return measurements::fixed_measurement(atan2(vel_[1].value(), vel_[0].value()), units::rad); 
                }     

                inline measurements::fixed_measurement theta() const { 
                    if (vel_[2] == 0) return measurements::fixed_measurement(0, units::rad);
                    return measurements::fixed_measurement(acos(vel_[2].value() / magnitude().value()), units::rad); 
                }

                inline std::vector<double> direction() const {
                    return { cos(phi().value()), sin(phi().value()), vel_[2].value() / magnitude().value() };
                } 

                void print_velocity() const {
                    std::cout << "velocity = {\n";
                    for (auto i : vel_) i.print(); 
                    std::cout << "}\n"; 
                }


        }; // class velocity
        

        namespace objects {


            class material_point : public position, public velocity {

                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================

                    explicit material_point(const position& pos = position(), const velocity& vel = velocity()) noexcept : 
                        position(pos), velocity(vel) {}

                    explicit material_point(const std::pair<position, velocity>& pos_vel) noexcept : 
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

                    void set_pos_vel(const std::pair<std::vector<measurements::fixed_measurement>, std::vector<measurements::fixed_measurement>>& pos_vel) {
                        set_position(pos_vel.first); 
                        set_velocity(pos_vel.second); 
                    }               

                    inline std::pair<std::vector<measurements::fixed_measurement>, 
                                    std::vector<measurements::fixed_measurement>> get_pos_vel() const { return std::make_pair(get_position(), get_velocity()); }

                    void print_pos_vel() const {
                        print_position(); 
                        print_velocity(); 
                    }


            }; // class material point


            class mass : virtual public material_point {

                protected: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::fixed_measurement mass_; 
                    

                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    explicit mass(const double& mass, const units::unit& unit = units::kg, const position& pos = position(), const velocity& vel = velocity()) noexcept : 
                        material_point(pos, vel), mass_{mass, unit} {}

                    explicit mass(const double& mass, const units::unit& unit, const std::pair<position, velocity>& pos_vel) noexcept :
                        material_point(pos_vel), mass_{mass, unit} {}

                    explicit mass(const measurements::fixed_measurement& mass, const position& pos = position(), const velocity& vel = velocity()) noexcept :
                        material_point(pos, vel), mass_{mass} {}

                    explicit mass(const measurements::fixed_measurement& mass, const std::pair<position, velocity>& pos_vel) noexcept :
                        material_point(pos_vel), mass_{mass} {}

                    mass(const mass& mass) noexcept : material_point(mass.get_pos_vel()), mass_{mass.mass_} {}

                    ~mass() = default; 


                    // =============================================
                    // operators
                    // =============================================

                    inline mass operator+(const mass& other) const { return mass(mass_ + other.mass_, get_pos_vel()); }                    
                    
                    inline mass operator-(const mass& other) const { return mass(mass_ - other.mass_, get_pos_vel()); }                    

                    inline mass operator+=(const mass& other) const { return mass(mass_ + other.mass_, get_pos_vel()); }                    
                    
                    inline mass operator-=(const mass& other) const { return mass(mass_ - other.mass_, get_pos_vel()); }                    
                    
                    inline mass operator*(const mass& other) const { return mass(mass_ * other.mass_, get_pos_vel()); }                    
                    
                    inline mass operator*(const double& val) const { return mass(mass_ * val, get_pos_vel()); }  

                    inline mass operator/(const mass& other) const { return mass(mass_ / other.mass_, get_pos_vel()); }                    

                    inline mass operator/(const double& val) const { return mass(mass_ / val, get_pos_vel()); }  


                    // =============================================
                    // set, get and print methods
                    // =============================================

                    constexpr void set_mass(const double& mass) { mass_.value(mass); }
                    
                    constexpr measurements::fixed_measurement get_mass() const { return mass_; }

                    void print_mass() const { 
                        std::cout << "mass = "; 
                        mass_.print(); 
                    }
                

            }; // class mass


            class charge : public virtual material_point {

                protected: 

                    // =============================================
                    // class member
                    // =============================================

                    measurements::fixed_measurement charge_; 
                    

                public: 

                    // =============================================
                    // constructors and destructor
                    // =============================================
                    
                    explicit charge(const double& charge, const units::unit& unit = units::kg, const position& pos = position(), const velocity& vel = velocity()) noexcept : 
                        material_point(pos, vel), charge_{charge, unit} {}

                    explicit charge(const double& charge, const units::unit& unit, const std::pair<position, velocity>& pos_vel) noexcept :
                        material_point(pos_vel), charge_{charge, unit} {}

                    explicit charge(const measurements::fixed_measurement& charge, const position& pos = position(), const velocity& vel = velocity()) noexcept :
                        material_point(pos, vel), charge_{charge} {}

                    explicit charge(const measurements::fixed_measurement& charge, const std::pair<position, velocity>& pos_vel) noexcept :
                        material_point(pos_vel), charge_{charge} {}

                    charge(const charge& charge) noexcept : material_point(charge.get_pos_vel()), charge_{charge.charge_} {}

                    ~charge() = default; 


                    // =============================================
                    // operators
                    // =============================================

                    inline charge operator+(const charge& other) const { return charge(charge_ + other.charge_, get_pos_vel()); }                    
                    
                    inline charge operator-(const charge& other) const { return charge(charge_ - other.charge_, get_pos_vel()); }                    

                    inline charge operator+=(const charge& other) const { return charge(charge_ + other.charge_, get_pos_vel()); }                    
                    
                    inline charge operator-=(const charge& other) const { return charge(charge_ - other.charge_, get_pos_vel()); }                    
                    
                    inline charge operator*(const charge& other) const { return charge(charge_ * other.charge_, get_pos_vel()); }                    
                    
                    inline charge operator*(const double& val) const { return charge(charge_ * val, get_pos_vel()); }  

                    inline charge operator/(const charge& other) const { return charge(charge_ / other.charge_, get_pos_vel()); }                    

                    inline charge operator/(const double& val) const { return charge(charge_ / val, get_pos_vel()); }  


                    // =============================================
                    // set & get & print methods
                    // =============================================

                    constexpr void set_charge(const double& charge) { charge_.value(charge); }
                    
                    constexpr measurements::fixed_measurement get_charge() const { return charge_; }

                    void print_charge() const { 
                        std::cout << "charge = "; 
                        charge_.print(); 
                    }

                
            }; // class charge


            class particle : public mass, public charge {

                public:  

                    // =============================================
                    // constructors & destructor
                    // =============================================     

                    explicit particle(const measurements::fixed_measurement& mass, const measurements::fixed_measurement& charge, const std::pair<position, velocity>& pos_vel) noexcept : 
                        material_point(pos_vel), mass::mass(mass), charge::charge(charge) {}

                    explicit particle(const measurements::fixed_measurement& mass, const measurements::fixed_measurement& charge, const position& pos = position(), const velocity& vel = velocity()) noexcept : 
                        material_point(pos, vel), mass::mass(mass), charge::charge(charge) {}

                    explicit particle(const mass& mass, const charge& charge, const material_point& pos_vel) noexcept : 
                        material_point(pos_vel), mass::mass(mass), charge::charge(charge) {}

                    explicit particle(const mass& mass, const charge& charge, const position& pos = position(), const velocity& vel = velocity()) noexcept : 
                        material_point(pos, vel), mass::mass(mass), charge::charge(charge) {}

                    particle(const particle& part) noexcept : particle(part.get_mass(), part.get_charge(), part.material_point::get_pos_vel()) {}

                    ~particle() = default; 


                    // =============================================
                    // print methods
                    // =============================================

                    void print_particle() const {
                        print_mass();
                        print_charge(); 
                        material_point::print_pos_vel();
                    }


            }; // class particle

        
        } // namespace objects


        // namespace defining some usefull constants
        namespace constants {

            objects::mass e_mass = objects::mass(9.109383701528e-31, units::kg);     
            objects::mass p_mass = objects::mass(1.672621923695e-27, units::kg);
            objects::charge e_charge = objects::charge(-1.602176634e-19, units::SI_derived::C);
            objects::charge p_charge = objects::charge(1.602176634e-19, units::SI_derived::C); 

        } // namespace constants


        namespace objects {

            particle electon = particle(constants::e_mass, constants::e_charge); 
            
            particle proton = particle(constants::p_mass, constants::p_charge);


            class time {

                public:     

                    // =============================================
                    // class members
                    // =============================================     

                    measurements::fixed_measurement time_; 


                    // =============================================
                    // constructor and destructor
                    // =============================================   

                    explicit constexpr time(const double& time = 0.0, const units::unit& unit = units::s) noexcept : time_(time, unit) {}
                    
                    ~time() = default; 

                    
                    // =============================================
                    // set & get & print methods
                    // =============================================   

                    constexpr void set_time(const double& t) { time_.value(t); }
                                        
                    constexpr measurements::fixed_measurement get_time() const { return time_; }

                    constexpr void increase_time(const double& t) { time_ += t; } 

                    constexpr void reset_time() { time_.value(0.0); }

                    inline void print_time() const { time_.print(); }


            }; // class time


        } // namespace objects


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

                timer(const units::unit& unit = units::s) noexcept : unit_{unit} {}

                ~timer() = default;

                
                // =============================================
                // timer methods
                // =============================================   

                inline void start() { start_ = std::chrono::system_clock::now(); }

                inline void pause() { pause_ = std::chrono::system_clock::now(); }
                
                void print() const { 
                    std::cout << "elapsed time = " << static_cast<std::chrono::duration<double>>(pause_ - start_).count() << " "; 
                    unit_.print(); 
                    std::cout << "\n";
                }


        }; // class timer


    } // namespace physics


    namespace math {


        namespace tools {


            class ode_solver : public physics::objects::time {

                    protected: 

                        // =============================================
                        // class members
                        // =============================================

                        std::pair<std::vector<physics::measurements::fixed_measurement>,
                                std::vector<physics::measurements::fixed_measurement>> m_df; 
                        
                        double h_{}; 


                    public: 

                        // =============================================
                        // virtual destructor
                        // =============================================

                        virtual ~ode_solver() {}


                        // =============================================
                        // virtual eval methods
                        // =============================================

                        virtual std::pair<std::vector<physics::measurements::fixed_measurement>,
                                        std::vector<physics::measurements::fixed_measurement>> eval(const std::pair<std::vector<physics::measurements::fixed_measurement>, std::vector<physics::measurements::fixed_measurement>>& pos_vel, const double& h = 0.001) = 0; 


                        // =============================================
                        // integration methods
                        // =============================================

                        inline std::pair<std::vector<physics::measurements::fixed_measurement>,
                                        std::vector<physics::measurements::fixed_measurement>> euler(const std::pair<std::vector<physics::measurements::fixed_measurement>, std::vector<physics::measurements::fixed_measurement>>& pos_vel, const double& h = 0.001) {
                            increase_time(h); 
                            return pos_vel + h * eval(pos_vel, h);
                        }

                        std::pair<std::vector<physics::measurements::fixed_measurement>,
                                std::vector<physics::measurements::fixed_measurement>> euler_modified(const std::pair<std::vector<physics::measurements::fixed_measurement>, std::vector<physics::measurements::fixed_measurement>>& pos_vel, const double& h = 0.001) {
                            std::pair<std::vector<physics::measurements::fixed_measurement>,
                                    std::vector<physics::measurements::fixed_measurement>> appo = pos_vel + h * eval(pos_vel, h); 
                            increase_time(h); 
                            return pos_vel + h * (appo + eval(pos_vel + h * appo, h)) / 2.; 
                        }

                        std::pair<std::vector<physics::measurements::fixed_measurement>,
                                std::vector<physics::measurements::fixed_measurement>> rk4(const std::pair<std::vector<physics::measurements::fixed_measurement>, std::vector<physics::measurements::fixed_measurement>>& pos_vel, const double& h = 0.001) {
                            std::pair<std::vector<physics::measurements::fixed_measurement>,
                                    std::vector<physics::measurements::fixed_measurement>> k1{}, k2{}, k3{}, k4{}; 
                            k1 = eval(pos_vel, h); 
                            k2 = eval(pos_vel + k1 * h / 2., h);
                            k3 = eval(pos_vel + k2 * h / 2., h);
                            k4 = eval(pos_vel + k3 * h / 2., h);      
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

                        measurements::fixed_measurement omega_;


                    public: 

                        // =============================================
                        // constructors and destructor
                        // =============================================
                        
                        harmonic(const measurements::fixed_measurement& omega, const measurements::fixed_measurement& time = measurements::fixed_measurement(0, units::s)) : omega_{omega} { time::time_ = time; }
                        
                        ~harmonic() {}


                        // =============================================
                        // set and get methods
                        // =============================================

                        constexpr void set_omega(const double& omega) { omega_ = omega; }

                        constexpr measurements::fixed_measurement get_omega() const { return omega_; }

                        void print_omega() const {
                            std::cout << "omega = "; 
                            omega_.print();
                        }


                        // =============================================
                        // eval methods
                        // =============================================

                        inline std::pair<std::vector<measurements::fixed_measurement>,
                                        std::vector<measurements::fixed_measurement>> eval(const std::pair<std::vector<physics::measurements::fixed_measurement>, std::vector<physics::measurements::fixed_measurement>>& pos_vel, const double& h = 0.001) override {
                            return std::make_pair(pos_vel.second, - math::algebra::square(omega_.value()) * pos_vel.first);
                        }


                }; // class harmonic


            } // namespace oscillators


        } // namespace objects


    } // namespace physics


} // namespace physim
