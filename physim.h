

// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physim(namespace) containing the basic tools for computational physics. 
// last updated:    20/08/2022


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
#include <unordered_map>
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

        // check if the two value given differ from each other less than the precision given
        constexpr bool are_close(const double& calculated, const double& expected, const double& epsilon = 1e-6){
            return std::fabs(calculated - expected) <= epsilon;
        }

        // class for keeping track of the inexorable passage of time
        class timer {

            public:

                // =============================================
                // class members
                // =============================================     
                
                std::chrono::duration<double> m_elapsed_seconds;
                std::chrono::time_point<std::chrono::system_clock> m_start, m_pause;


                // =============================================
                // constructor and destructor
                // =============================================   

                // timer(const units::unit& unit) : time(0.0, unit) {}

                // timer() : time(0.0, units::defined::s) {}

                timer() {}

                ~timer() = default;

                
                // =============================================
                // timer methods
                // =============================================   

                void start() { m_start = std::chrono::system_clock::now(); }

                void pause() { m_pause = std::chrono::system_clock::now(); }

                double elapsed_time() { 
                    m_elapsed_seconds = m_pause - m_start;
                    return m_elapsed_seconds.count(); 
                }
                
                void print() { std::cout << "elapsed time = " << elapsed_time() << "\n"; }

        }; // class timer


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
                        if (value < T{0.0}) { return constants::invalid_conversion; }
                        else return std::sqrt(value);
                    case -2:
                        if (value < T{0.0}) { return constants::invalid_conversion; }
                        else return std::sqrt(T{1.0} / value);
                    case 3: 
                        return std::cbrt(value);
                    case -3: 
                        return std::cbrt(T{1.0} / value);
                    case 4: 
                        if (value < T{0.0}) { return constants::invalid_conversion; }
                        else return std::sqrt(std::sqrt(value));
                    case -4: 
                        if (value < T{0.0}) { return constants::invalid_conversion; }
                        else return std::sqrt(std::sqrt(T{1.0} / value));
                    default:
                        if (value < T{0.0} && power % 2 == 0) { return constants::invalid_conversion; }
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
                        std::cout << "f (" << x << ") = " << std::setprecision((int)-log10(precision)) << eval(x) << "\n"; 
                    }
                
                    virtual void print_equation() const = 0;
                

            }; // class function_base
            
            // functor class for composing single functions into some more complex expressions
            class functor : public function_base {

                private:

                    // =============================================
                    // class members
                    // =============================================
                
                    char m_op; 
                    function_base * m_f; 
                    function_base * m_g;
                
                
                public:

                    // =============================================
                    // constructor and destructor
                    // =============================================     
                
                    functor(const char& op, function_base* f, function_base* g) {
                        m_op = op;
                        m_f = f; 
                        m_g = g;
                    }
                    
                    ~functor() {}
                
                
                    // =============================================
                    // eval methods
                    // =============================================
                
                    double eval(const double& x) const override {
                        switch (m_op) {
                            case '+':
                                return m_f->eval(x) + m_g->eval(x);

                            case '-': 
                                return m_f->eval(x) - m_g->eval(x);

                            case '*':
                                return m_f->eval(x) * m_g->eval(x);

                            case '/': 
                                return m_f->eval(x) / m_g->eval(x);
                            
                            case '^':
                                return std::pow(m_f->eval(x), m_g->eval(x));

                            case 'c':
                                return m_f->eval(m_g->eval(x));
                            
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
                    double m_m, m_q; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================     
            
                    explicit line(const double& m = 1, const double& q = 0) noexcept : m_m{m}, m_q{q} {}
                    
                    ~line() {}

            
                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void m(const double& m) { m_m = m; } 

                    constexpr void q(const double& q) { m_q = q; }
            
                    constexpr double m() const { return m_m; }

                    constexpr double q() const { return m_q; } 


                    // =============================================
                    // eval methods
                    // =============================================

                    constexpr double eval(const double& x) const override { return m_m * x + m_q; }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_m != 1 && m_m != -1) std::cout << m_m; 
                        else {
                            if (m_m == 1) std::cout << "+";
                            if (m_m == -1) std::cout << "-";
                        }
                        std::cout << "x"; 
                        if (m_q != 0) {
                            if (m_q > 0) std::cout << " + " << m_q << "\n"; 
                            else std::cout << " - " << std::fabs(m_q) << "\n"; 
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
                    double m_a, m_b, m_c, m_delta; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================     
            
                    explicit quadratic(const double& a, const double& b, const double& c) noexcept 
                        : m_a{a}, m_b{b}, m_c{c}, m_delta{math::algebra::square(m_b) - 4 * m_a * m_c} {}
                    
                    ~quadratic() {}

            
                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void a(const double& a) { m_a = a; } 

                    constexpr void b(const double& b) { m_b = b; }

                    constexpr void c(const double& c) { m_c = c; }  
            
                    constexpr double a() const { return m_a; }

                    constexpr double b() const { return m_b; } 

                    constexpr double c() const { return m_c; } 

                    constexpr double delta() const { return m_delta; } 

                    std::pair<double, double> roots() const {
                        if (m_delta == 0) return std::make_pair(- m_b / (2 * m_a), - m_b / (2 * m_a));
                        else if (m_delta > 0) return std::make_pair(- m_b - algebra::sqrt(m_delta) / (2 * m_a), (- m_b + algebra::sqrt(m_delta)) / (2 * m_a));
                        else std::cout << "\nThere are not real solutions...\n\n"; // note: create a complex class then come back here
                        return std::make_pair(NAN, NAN); 
                    }


                    // =============================================
                    // eval methods
                    // =============================================

                    constexpr double eval_Horner(const double& x) const { return m_c + x * (m_b + x * m_a); }

                    constexpr double eval(const double& x) const override { return m_a * algebra::square(x) + m_b * x + m_c; }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_a != 1 && m_a != -1) std::cout << m_a; 
                        else if (m_a == -1) std::cout << "-"; 
                        else std::cout << "x^2 "; 
                        if (m_b != 0) {
                            if (m_b > 0) std::cout << "+ ";
                            else std::cout << "- ";
                            if (m_b != 1 && m_b != -1) std::cout << std::fabs(m_b) << "x ";  
                            else std::cout << "x "; 
                        } 
                        if (m_c != 0) {
                            if (m_c > 0) std::cout << "+ ";
                            else std::cout << "- "; 
                            std::cout << std::fabs(m_c) << "\n";
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
                    double m_a, m_b, m_c, m_d; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
                
                    explicit cubic(const double& a, const double& b, const double& c, const double& d) noexcept
                        : m_a{a}, m_b{b}, m_c{c}, m_d{d} {}
                
                    ~cubic() {}  
                
                
                    // =============================================
                    // set methods
                    // =============================================
            
                    constexpr void a(const double& a) { m_a = a; } 

                    constexpr void b(const double& b) { m_b = b; }

                    constexpr void c(const double& c) { m_c = c; }  
            
                    constexpr void d(const double& d) { m_d = d; } 
            
                
                    // =============================================
                    // get methods
                    // =============================================

                    constexpr double a() const { return m_a; }

                    constexpr double b() const { return m_b; } 

                    constexpr double c() const { return m_c; } 

                    constexpr double d() const { return m_d; } 

                
                    // =============================================
                    // eval methods
                    // =============================================

                    constexpr double eval(const double& x) const override { 
                        return m_a * math::algebra::cube(x) + m_b * math::algebra::square(x) + m_c * x + m_d; 
                    }
            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_a != 1) std::cout << m_a; 
                        std::cout << "x^3 "; 
                        if (m_b != 0) {
                            if (m_b != 1) std::cout << "+ " << m_b << "x^2 "; 
                            else std::cout << "+ x^2 "; 
                        }
                        if (m_c != 0) {
                            if (m_c != 1) std::cout << "+ " << m_c << "x"; 
                            else std::cout << "+ x"; 
                        }
                        if (m_d != 0) std::cout << " + " << m_d << "\n";
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
                    double m_c; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
                
                    explicit square_root(const double& c = 1) noexcept : m_c{c} {}
                
                    ~square_root() {}
                
                
                    // =============================================
                    // set & get methods
                    // =============================================
                
                    constexpr void c(const double& c) { m_c = c; }

                    constexpr double c() const { return m_c; }

                
                    // =============================================
                    // eval methods
                    // =============================================

                    double eval(const double& x) const override { return m_c * math::algebra::sqrt(x); }
                            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c != 1) std::cout << m_c; 
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
                    double m_c; 


                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
                
                    explicit cubic_root(double c = 1) noexcept : m_c{c} {}
                
                    ~cubic_root() {} 
                
                
                    // =============================================
                    // set & get methods
                    // =============================================
                
                    constexpr void c(double c) { m_c = c; }

                    constexpr double c() const { return m_c; }

                
                    // =============================================
                    // eval methods
                    // =============================================

                    double eval(const double& x) const override { return m_c * math::algebra::cbrt(x); }
                            
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c != 1) std::cout << m_c; 
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
                    double m_base, m_c1, m_c2;

                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
                    explicit exponential(const double& base = math::constants::e, const double& c1 = 1, const double& c2 = 1) noexcept
                        : m_base{base}, m_c1{c1}, m_c2{c2} {}

                    ~exponential() {} 
            

                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void base(const double& base) { m_base = base; } 
                
                    constexpr void c1(const double& c1) { m_c1 = c1; }
                
                    constexpr void c2(const double& c2) { m_c2 = c2; }
            
                    constexpr double base() const { return m_base; }
            
                    constexpr double c1() const { return m_c1; }

                    constexpr double c2() const { return m_c2; }
                
                
                    // =============================================
                    // eval methods
                    // =============================================
            
                    constexpr double eval(const double& x) const override { return m_c1 * std::pow(m_base, m_c2 * x); }
                            
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        if (m_base != math::constants::e) std::cout << m_base << "^"; 
                        else std::cout << "e^(";
                        if (m_c2 != 1 && m_c2 != -1) std::cout << m_c2; 
                        if (m_c2 == -1) std::cout << "-"; 
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
                    double m_base, m_c1, m_c2; 

                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
                    explicit logarithm(const double& base = math::constants::e, const double& c1 = 1, const double& c2 = 1) noexcept
                        : m_base{base}, m_c1{c1}, m_c2{c2} {}
                
                    ~logarithm() {} 
            

                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void base(const double& base) { m_base = base; } 
            
                    constexpr void c1(const double& c1) { m_c1 = c1; }
                
                    constexpr void c2(const double& c2) { m_c2 = c2; }
                            
                    constexpr double base() const { return m_base; }

                    constexpr double c1() const { return m_c1; }

                    constexpr double c2() const { return m_c2; }
                

                    // =============================================
                    // eval methods
                    // =============================================
            
                    constexpr double eval(const double& x) const override { return m_c1 * std::log(m_c2 * x) / std::log(m_base); }
                                
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "log";
                        if (m_base != math::constants::e) std::cout << "_( " << m_base << ")";
                        if (m_c2 != 1) std::cout << "(" << m_c2 << "x) \n"; 
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
                    double m_c1, m_c2; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
                    explicit sine(const double& c1 = 1, const double& c2 = 1) noexcept : m_c1{c1}, m_c2{c2} {} 
                
                    ~sine() {} 

                
                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void c1(const double& c1) { m_c1 = c1; }
                
                    constexpr void c2(const double& c2) { m_c2 = c2; }
                
                    constexpr double c1() const { return m_c1; }

                    constexpr double c2() const { return m_c2; }
                
                
                    // =============================================
                    // eval methods
                    // =============================================
            
                    constexpr double eval(const double& x) const override { return m_c1 * std::sin(m_c2 * x); }
                        
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "sin(";
                        if (m_c2 != 1) std::cout << m_c2; 
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
                    double m_c1, m_c2; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
            
                    explicit cosine(const double& c1 = 1, const double& c2 = 1) noexcept : m_c1{c1}, m_c2{c2} {} 
                
                    ~cosine() {} 

                
                    // =============================================
                    // set & get methods
                    // =============================================
            
                    constexpr void c1(const double& c1) { m_c1 = c1; }
                
                    constexpr void c2(const double& c2) { m_c2 = c2; }
                
                    constexpr double c1() const { return m_c1; }

                    constexpr double c2() const { return m_c2; }
                

                    // =============================================
                    // eval methods
                    // =============================================
            
                    constexpr double eval(const double& x) const override { return m_c1 * std::cos(m_c2 * x); }
                        
            
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "cos(";
                        if (m_c2 != 1) std::cout << m_c2; 
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
                    double m_c1, m_c2; 
                
                
                public: 
            
                    // =============================================
                    // constructor and destructor
                    // =============================================   
            
                    explicit tangent(const double& c1 = 1, const double& c2 = 1) noexcept : m_c1{c1}, m_c2{c2} {} 
                
                    ~tangent() {} 

                
                    // =============================================
                    // set methods
                    // =============================================
            
                    constexpr void c1(const double& c1) { m_c1 = c1; }
                
                    constexpr void c2(const double& c2) { m_c2 = c2; }
                
                
                    // =============================================
                    // get methods
                    // =============================================

                    constexpr double c1() const { return m_c1; }

                    constexpr double c2() const { return m_c2; }
                

                    // =============================================
                    // eval methods
                    // =============================================
            
                    constexpr double eval(const double& x) const override { return m_c1 * std::tan(m_c2 * x); }
                    
            
                    // =============================================
                    // print methods
                    // =============================================

                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        std::cout << "tan(";
                        if (m_c2 != 1) std::cout << m_c2; 
                        std::cout << "x) \n";
                    }


            }; // class tangent


        } // namespace functions
        

        // class random_generator for generating pseudo-casual numbers
        class random_generator {

            private: 

                // =============================================
                // class members
                // =============================================        
                
                unsigned int m_a, m_c, m_m, m_seed;


            public: 
        
                // =============================================
                // constructors
                // =============================================
                
                random_generator() : 
                    m_a{1664525}, m_c{1013904223}, m_m{static_cast<unsigned int>(std::pow(2, 31))} {}

                random_generator(const unsigned int& seed) 
                    : m_a{1664525}, m_c{1013904223}, m_m{static_cast<unsigned int>(std::pow(2, 31))}, m_seed{seed} {}


                // =============================================
                // set & get methods
                // =============================================

                constexpr void a(const unsigned int& a) { m_a = a; }
                
                constexpr void c(const unsigned int& c) { m_c = c; }
                
                constexpr void m(const unsigned int& m) { m_m = m; }
                
                constexpr void seed(const unsigned int& seed) { m_seed = seed; }

                constexpr unsigned int a() const { return m_a; }
                
                constexpr unsigned int c() const { return m_c; }
                
                constexpr unsigned int m() const { return m_m; }
                
                constexpr unsigned int seed() const { return m_seed; }


                // =============================================
                // distributions methods
                // =============================================

                constexpr double rand(const double& min = 0., const double& max = 1.) {
                    seed(static_cast<unsigned int>((m_a * m_seed + m_c) % m_m)); 
                    return min + std::fabs(max - min) * m_seed / m_m; 
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

        }; // class random_generator


        // class integral for evaluating the integrals
        class integral {

            private: 

                // =============================================
                // class members
                // =============================================

                double m_a{}, m_b{}, m_h{};
                
                unsigned int m_steps{};

                int m_sign{}; 

                double m_sum{}, m_integral{}, m_old_integral{}, m_error{};  

                math::random_generator m_rg;


                // =============================================
                // set methods
                // =============================================

                constexpr void steps(const unsigned int& n) { 
                    m_steps = n; 
                    m_h = std::fabs(m_b - m_a) / m_steps;
                }        

                constexpr void check_range() { m_sign = (m_a == m_b ? 0. : (m_b > m_a ? 1. : -1)); }
                
                constexpr void sum(const double& sum) { m_sum = sum; }

                constexpr void reset_integral() { m_integral = 0; }    

                constexpr void begin_integration(const double& a, const double& b, unsigned int n = 1000, const double& sum0 = 0) {
                    m_a = a; 
                    m_b = b; 
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

                constexpr double a() const { return m_a; }

                constexpr double b() const { return m_b; }

                constexpr int sign() const { return m_sign; }

                constexpr unsigned int steps() const { return m_steps; }

                constexpr double h() const { return m_h; }

                constexpr double sum() const { return m_sum; }

                constexpr double value() const { return m_integral; }

                constexpr double error() const { return m_error; }


                // =============================================
                // print methods
                // =============================================

                constexpr void error_integral(const double& delta) { m_error = 4 * delta / 3.; } 

                void print_value(const double& precision = 1.e-6) {
                    std::cout << "\nIntegral of f(x) in [" << m_a << ", " << m_b << "] = " << std::setprecision((int)-log10(precision)) << m_integral << "\n";
                }

                void print_error(const double& precision = 1.e-6) {
                    std::cout << "error = " << std::setprecision((int)-log10(precision)) << m_error << "\n";
                }        

                void print(const double& precision = 1.e-6) {
                    print_value(); 
                    print_error(); 
                }
                
                
                // =============================================
                // integration methods
                // =============================================

                void midpoint(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    for (unsigned int i{}; i < m_steps; i++) { m_sum += (f.eval(m_a + (i + 0.5) * m_h)); }
                    m_integral = m_sum * m_h; 
                }

                void midpoint_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    begin_integration(a, b, 1); 
                    while (true) {
                        m_old_integral = m_integral; 
                        midpoint(m_a, m_b, f, m_steps * 2);
                        error_integral(std::fabs(m_integral - m_old_integral));
                        if (m_error < prec) break;
                    }    
                }
                
                void trapexoid(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 2.);
                    for (unsigned int i{1}; i < m_steps; i++) m_sum += f.eval(m_a + i * m_h); 
                    m_integral = m_sum * m_h; 
                } 

                void trapexoid_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    begin_integration(a, b, 2, f.eval(a) + f.eval(b) / 2. + f.eval((a + b) / 2.)); 
                    while (true) {
                        m_old_integral = m_integral; 
                        trapexoid(m_a, m_b, f, m_steps * 2);
                        error_integral(std::fabs(m_integral - m_old_integral));
                        if (m_error < prec) break;
                    }
                }

                void simpson(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    if (n % 2 == 0) begin_integration(a, b, n, (f.eval(a) + f.eval(b)) / 3.);
                    else begin_integration(a, b, n + 1);  
                    for (unsigned int i{1}; i < m_steps; i++) m_sum += 2 * (1 + i % 2) * (f.eval(m_a + i * m_h)) / 3.; 
                    m_integral = m_sum * m_h; 
                }

                void simpson_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    begin_integration(a, b, 2, (f.eval(a) + f.eval(b)) / 3.); 
                    while (true) {
                        m_old_integral = m_integral; 
                        simpson(m_a, m_b, f, m_steps * 2);
                        error_integral(std::fabs(m_integral - m_old_integral));
                        if (m_error < prec) break; 
                    }
                }

                void mean(const double& a, const double& b, const functions::function_base& f, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    for (unsigned int i{}; i < n; i ++) m_sum += f.eval(m_rg.rand(a, b)); 
                    m_integral = (m_b - m_a) * m_sum / m_steps; 
                }

                void mean_fixed(const double& a, const double& b, const functions::function_base& f, const double& prec = 1.e-6) {
                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        mean(a, b, f);
                        k.emplace_back(m_integral); 
                    }
                    double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                    mean(a, b, f, static_cast<unsigned int>(math::algebra::square(k_mean / prec))); 
                }
        
                void hit_or_miss(const double& a, const double& b, const functions::function_base& f, const double& fmax, unsigned int n = 1000) {
                    begin_integration(a, b, n); 
                    double x{}, y{}; 
                    unsigned int hits{};
                    for (unsigned int i{}; i < n; i ++) {
                        x = m_rg.rand(a, b); 
                        y = m_rg.rand(0., fmax);  
                        if (y <= f.eval(x)) hits++; 
                    }
                    m_integral = (m_b - m_a) * fmax * hits / n; 
                }

                void hit_or_miss_fixed(const double& a, const double& b, const functions::function_base& f, const double& fmax, const double& prec = 1.e-6) {
                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        hit_or_miss(a, b, f, fmax);
                        k.emplace_back(m_integral); 
                    }
                    double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                    unsigned int N = static_cast<unsigned int>(math::algebra::square(k_mean / prec)); 
                    hit_or_miss(a, b, f, fmax, N); 
                }


        }; // class integral


    } // namespace math


} // namespace physim
