

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
            
                    explicit quadratic(const double& a, const double& b, const double& c) noexcept :
                        m_a{a}, m_b{b}, m_c{c}, m_delta{math::algebra::square(m_b) - 4 * m_a * m_c} {}
                    
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
                
                    explicit cubic(const double& a, const double& b, const double& c, const double& d) noexcept :
                        m_a{a}, m_b{b}, m_c{c}, m_d{d} {}
                
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
            
                    explicit exponential(const double& base = math::constants::e, const double& c1 = 1, const double& c2 = 1) noexcept :
                        m_base{base}, m_c1{c1}, m_c2{c2} {}

                    explicit exponential(const double& c1 = 1, const double& c2 = 1) noexcept :
                        m_base{math::constants::e}, m_c1{c1}, m_c2{c2} {}

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
            
                    explicit logarithm(const double& base = math::constants::e, const double& c1 = 1, const double& c2 = 1) noexcept : 
                        m_base{base}, m_c1{c1}, m_c2{c2} {}

                    explicit logarithm(const double& c1 = 1, const double& c2 = 1) noexcept : 
                        m_base{math::constants::e}, m_c1{c1}, m_c2{c2} {}

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
        

        // namespace defining some usefull operations
        namespace tools {
            
            // check if the two value given differ from each other less than the precision given
            constexpr bool are_close(const double& calculated, const double& expected, const double& epsilon = 1e-6){
                return std::fabs(calculated - expected) <= epsilon;
            }

            // round a value to the expected level of precision of a double
            inline double cround(double val) {
                std::uint64_t bits;
                std::memcpy(&bits, &val, sizeof(bits));
                bits += 0x800ULL;
                bits &= 0xFFFFFFFFFFFFF000ULL;
                std::memcpy(&val, &bits, sizeof(bits));
                return val;
            }

            // rounding compare for equality on double
            inline bool compare_round_equals(double val1, double val2) {
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
                    
                    unsigned int m_a, m_c, m_m, m_seed;


                public: 
            
                    // =============================================
                    // constructors & destructor
                    // =============================================
                    
                    random_generator() : 
                        m_a{1664525}, m_c{1013904223}, m_m{static_cast<unsigned int>(std::pow(2, 31))} {}

                    random_generator(const unsigned int& seed) :
                        m_a{1664525}, m_c{1013904223}, m_m{static_cast<unsigned int>(std::pow(2, 31))}, m_seed{seed} {}

                    ~random_generator() {}

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

            }; // class tools::random_generator


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

                    random_generator m_rg;


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
      
                    const char* name_ = "-";  


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

                    static constexpr uint32_t bits[7] = { 
                        bitwidth::metre, bitwidth::second, 
                        bitwidth::kilogram, bitwidth::ampere, 
                        bitwidth::kelvin, bitwidth::mole, 
                        bitwidth::candela
                    };  


                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // constructor from powers
                    constexpr unit_data(const int& metres, const int& seconds, const int& kilograms, const int& amperes, const int& kelvins,  const int& moles, const int& candelas) :
                        metre_(metres), second_(seconds), kilogram_(kilograms),
                        ampere_(amperes), kelvin_(kelvins), mole_(moles), candela_(candelas) {};

                    // constructor from powers and name
                    constexpr unit_data(const int& metres, const int& seconds, const int& kilograms, const int& amperes, const int& kelvins,  const int& moles, const int& candelas, const char* name_) :
                        metre_(metres), second_(seconds), kilogram_(kilograms), ampere_(amperes), 
                        kelvin_(kelvins), mole_(moles), candela_(candelas), name_(name_) {};

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
                    
                    // comparison operators
                    constexpr bool operator==(const unit_data& other) const {
                        return metre_ == other.metre_ && second_ == other.second_ &&
                            kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                            candela_ == other.candela_ && kelvin_ == other.kelvin_ && mole_ == other.mole_;         
                    }

                    // definitely not a comparison operators
                    constexpr bool operator!=(const unit_data& other) const {
                        return !(*this == other);
                    }

                    // check if the units have the same base unit 
                    constexpr bool has_same_base(const unit_data& other) const {
                        return *this == other;
                    }

                    // check if the unit is empty
                    constexpr bool empty() const {
                        return metre_ == 0 && second_ == 0 && kilogram_ == 0 &&
                            ampere_ == 0 && candela_ == 0 && kelvin_ == 0 && mole_ == 0;
                    }
                    

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

                    constexpr const char* name() const { return name_; }

                    constexpr int unit_type_count() const {
                        return ((metre_ != 0) ? 1 : 0) + ((second_ != 0) ? 1 : 0) +
                            ((kilogram_ != 0) ? 1 : 0) + ((ampere_ != 0) ? 1 : 0) +
                            ((kelvin_ != 0) ? 1 : 0) + ((mole_ != 0) ? 1 : 0) + ((candela_ != 0) ? 1 : 0);
                    }

                    constexpr void print() const {
                        if (name_[0] == '-') {
                            if (metre_ != 0 && metre_ != 1) std::cout << "m^" << metre_; 
                            if (metre_ == 1) std::cout << "m";
                            if (second_ != 0 && second_ != 1) std::cout << "s^" << second_; 
                            if (second_ == 1) std::cout << "s"; 
                            if (kilogram_ != 0 && kilogram_ != 1) std::cout << "kg^" << kilogram_; 
                            if (kilogram_ == 1) std::cout << "kg"; 
                            if (ampere_ != 0 && ampere_ != 1) std::cout << "A^" << ampere_; 
                            if (ampere_ == 1) std::cout << "A"; 
                            if (kelvin_ != 0 && kelvin_ != 1) std::cout << "K^" << kelvin_; 
                            if (kelvin_ == 1) std::cout << "K";
                            if (mole_ != 0 && mole_ != 1) std::cout << "mol^" << mole_; 
                            if (mole_ == 1) std::cout << "mol"; 
                            if (candela_ != 0 && candela_ != 1) std::cout << "cd^" << candela_; 
                            if (candela_ == 1) std::cout << "cd"; 
                        }
                        else std::cout << name_ << "\n"; 
                    }


                private: 

                    // check if the base_unit has a valid root
                    constexpr bool has_valid_root(const int& power) const {
                        return metre_ % power == 0 && second_ % power == 0 &&
                            kilogram_ % power == 0 && ampere_ % power == 0 &&
                            candela_ % power == 0 && kelvin_ % power == 0 &&
                            mole_ % power == 0;
                    }      
                    
                    // to handle a few weird operations that operate on square_root Hz
                    constexpr int root_Hertz_modifier(const int& power) const {
                        return (second_ * power == 0 || power % 2 != 0) ? 0 : (power / 2) * ((second_ < 0) || (power < 0) ? 9 : -9);
                    }    


            }; // class unit_base


            // class defining a basic unit module 
            class unit_prefix {

                public:  

                    // =============================================
                    // class members
                    // ============================================= 
                    
                    double multiplier_{1.0};  


                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    unit_prefix() noexcept {}
                    
                    explicit constexpr unit_prefix(const double& mult) noexcept : multiplier_{mult} {} 


                    // =============================================
                    // operators
                    // ============================================= 

                    // perform a multiply operation by adding the powers together
                    constexpr unit_prefix operator*(const unit_prefix& other) const {
                        return unit_prefix(multiplier_ * other.multiplier_);
                    }

                    // perform a division operation by subtract the powers together
                    constexpr unit_prefix operator/(const unit_prefix& other) const {
                        return unit_prefix(multiplier_ / other.multiplier_);
                    }

                    // invert the unit
                    constexpr unit_prefix inv() const {
                        return unit_prefix(1 / multiplier_);
                    }

                    // take a unit_prefix to some power
                    constexpr unit_prefix pow(const int& power) const { 
                        return unit_prefix(std::pow(multiplier_, power));
                    }
                    
                    // take some root of a unit_prefix
                    unit_prefix root(const int& power) const {
                        return unit_prefix(math::algebra::root(multiplier_, power)); 
                    }
                    
                    
                    // =============================================
                    // check methods
                    // ============================================= 

                    // comparison operators
                    constexpr bool operator==(const unit_prefix& other) const {
                        if (multiplier_ == other.multiplier_) return true;    
                        else return math::tools::compare_round_equals(multiplier_, other.multiplier_);
                    }

                    // definitely not a comparison operators
                    constexpr bool operator!=(const unit_prefix& other) const {
                        return !(*this == other);
                    }

                    // check if the units have the same prefix unit 
                    constexpr bool has_same_prefix(const unit_prefix& other) const {
                        return *this == other;
                    }

                    // check if the multiplier is nan
                    constexpr bool is_nan(const unit_prefix& pref) {
                        return std::isnan(pref.multiplier_);
                    }

                    // checks that the multiplier is finite
                    constexpr bool is_finite(const unit_prefix& pref) {
                        return std::isfinite(pref.multiplier_);
                    }

                    // check if the multiplier is infinite
                    constexpr bool is_inf(const unit_prefix& pref) {
                        return std::isinf(pref.multiplier_);
                    }


                    // =============================================
                    // get & print methods
                    // ============================================= 

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
            class unit {

                public:

                    // =============================================
                    // class members
                    // ============================================= 

                    unit_data data_{0, 0, 0, 0, 0, 0, 0}; 

                    unit_prefix prefix_{}; 


                    // =============================================
                    // constructors
                    // ============================================= 
                    
                    // default constructor
                    unit() noexcept {}

                    // constructor from unit_data 
                    explicit unit(const unit_data& unit_base) noexcept : data_{unit_base} {}

                    // constructor from unit_data and unit_prefix
                    explicit constexpr unit(const unit_data& unit_base, const unit_prefix& unit_prefix) noexcept : 
                        data_{unit_base}, prefix_{unit_prefix} {}

                    // constructor from unit_prefix and unit_data
                    explicit constexpr unit(const unit_prefix& unit_prefix, const unit_data& unit_base) noexcept : 
                        data_{unit_base}, prefix_{unit_prefix} {}

                    // constructor from unit_data and a double
                    constexpr unit(const unit_data& unit_base, const double& mult) noexcept : 
                        data_{unit_base}, prefix_{mult} {}

                    // constructor from double with a unit_data
                    constexpr unit(const double& mult, const unit_data& unit_base) noexcept :
                        data_{unit_base}, prefix_{mult} {}

                    // constructor from unit
                    constexpr unit(const unit& other) noexcept : 
                        data_{other.data_}, prefix_{other.prefix_} {} 

                    // constructor from unit and unit_prefix
                    constexpr unit(const unit& other, const unit_prefix& prefix) noexcept : 
                        data_{other.data_}, prefix_{prefix_} {} 

                    // constructor from unit_prefix and unit
                    constexpr unit(const unit_prefix& prefix, const unit& other) noexcept : 
                        data_{other.data_}, prefix_{prefix} {} 

                    // constructor from unit and a double
                    constexpr unit(const unit& other, const double& mult) noexcept : 
                        data_{other.data_}, prefix_{mult} {} 

                    // constructor from double with a unit
                    constexpr unit(const double& mult, const unit& other) noexcept :
                        data_{other.data_}, prefix_{mult} {} 


                    // =============================================
                    // operations
                    // ============================================= 

                    // take the reciprocal of a unit
                    constexpr unit inv() const { return unit(data_.inv(), prefix_.inv()); }

                    // multiply with another unit
                    constexpr unit operator*(const unit& other) const {
                        return unit(data_ * other.data_, prefix_ * other.prefix_);
                    }

                    // division operator
                    constexpr unit operator/(const unit& other) const {
                        return unit(data_ / other.data_, prefix_ / other.prefix_);
                    }                    

                    // take a unit to a power
                    constexpr unit pow(const int& power) const {
                        return unit(data_.pow(power), prefix_.pow(power));
                    }

                    // take a unit to a root power
                    unit root(const int& power) const {
                        return unit(data_.root(power), prefix_.root(power));                        
                    }


                    // =============================================
                    // checks
                    // ============================================= 

                    // equality operator
                    constexpr bool operator==(const unit& other) const {
                        if (data_ != other.data_) { return false; }
                        if (prefix_ == other.prefix_) { return true; }
                        return math::tools::compare_round_equals(multiplier(), other.multiplier());
                    }

                    // equality operator
                    constexpr bool operator!=(const unit& other) const { return !operator==(other); }

                    // test for exact numerical equivalence
                    constexpr bool is_exactly_the_same(const unit& other) const {
                        return data_ == other.data_ && prefix_ == other.prefix_;
                    }


                    // =============================================
                    // get & print methods
                    // ============================================= 

                    // get the multiplier
                    constexpr double multiplier() const { return prefix_.multiplier_; }
                    
                    // get the unit_data
                    constexpr unit_data base_units() const { return data_; }

                    // get the unit
                    constexpr unit as_unit() const { return *this; }

                    // print the unit
                    constexpr void print() const {
                        prefix_.print(); 
                        data_.print();
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
            namespace derived_units {

                // unit hertz(unit_data(0, 0, -1, 0, 0, 0, 0), "Hz");
                // unit Hz = hertz;
                // unit volt(unit_data(2, 1, -3, -1, 0, 0, 0), "V");
                // unit V = volt;
                // unit newton(unit_data(1, 1, -2, 0, 0, 0, 0), "N");
                // unit N = newton;
                // unit Pa(unit_data(-1, 1, -2, 0, 0, 0, 0), "Pa");
                // unit pascal = Pa;
                // unit joule(unit_data(2, 1, -2, 0, 0, 0, 0), "J");
                // unit J = joule;
                // unit watt(unit_data(2, 1, -3, 0, 0, 0, 0), "W");
                // unit W = watt;
                // unit coulomb(unit_data(0, 0, 1, 1, 0, 0, 0), "C");
                // unit C = coulomb;
                // unit farad(unit_data(-2, -1, 4, 2, 0, 0, 0), "F");
                // unit F = farad;
                // unit weber(unit_data(2, 1, -2, -1, 0, 0, 0), "Wb");
                // unit Wb = weber;
                // unit tesla(unit_data(0, 1, -2, -1, 0, 0, 0), "T");
                // unit T = tesla;
                // unit henry(unit_data(2, 1, -2, -2, 0, 0, 0), "H");                    
                // unit H = henry;
                unit mps(m / s);
                unit mpss(m / s.pow(2)); 

            } // namespace derived_units
            

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

            // generate a conversion factor between two physics::units::units in a constexpr function, the
            // physics::units::units will only convert if they have the same base physics::units::unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(UX start, UX2 result) { return quick_convert(1.0, start, result); }

            // generate a conversion factor between two physics::units::units in a constexpr function, the physics::units::units will only convert if they have the same base physics::units::unit
            template<typename UX, typename UX2>
            constexpr double quick_convert(double val, const UX& start, const UX2& result) {
                static_assert(std::is_same<UX, physics::units::unit>::value || std::is_same<UX, physics::units::unit>::value,
                    "convert argument types must be physics::units::unit or physics::units::unit");
                static_assert(std::is_same<UX2, physics::units::unit>::value || std::is_same<UX2, physics::units::unit>::value,
                    "convert argument types must be physics::units::unit or physics::units::unit");
                return (start.base_units() == result.base_units()) ?  val * start.multiplier() / result.multiplier() : 
                    constants::invalid_conversion;
            }

            // convert a value from one physics::units::unit base to another
            template<typename UX, typename UX2>
            double convert(double val, const UX& start, const UX2& result) {
                static_assert(std::is_same<UX, physics::units::unit>::value || std::is_same<UX, physics::units::unit>::value,
                    "convert argument types must be physics::units::unit or physics::units::unit");
                static_assert(std::is_same<UX2, physics::units::unit>::value || std::is_same<UX2, physics::units::unit>::value, 
                    "convert argument types must be physics::units::unit or physics::units::unit");     
                if (start == result) { return val; }
                if (start.base_units() == result.base_units()) { return val * start.multiplier() / result.multiplier(); }
                auto base_start = start.base_units();
                auto base_result = result.base_units();
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
                    measurement() noexcept {};

                    // constructor from a value and a unit 
                    explicit constexpr measurement(const double& val, const units::unit& unit) noexcept :
                        value_(val), unit_(unit) {}

                    // constructor from a measurement
                    constexpr measurement(const measurement& other) :
                        value_(other.value_), unit_(other.unit_) {}
                        
                    // destructor
                    ~measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    measurement& operator=(const measurement& val) noexcept {
                        value_ = (unit_ == val.units()) ? val.value() : val.value_as(unit_);
                        return *this;
                    }

                    measurement& operator=(const double& val) noexcept {
                        value_ = val;
                        return *this;
                    }     
    
                    constexpr measurement operator*(const measurement& other) const {
                        return measurement(value_ * other.value_, unit_ * other.unit_);
                    }
                
                    constexpr measurement operator*(const double& val) const {
                        return measurement(value_ * val, unit_);
                    }
                    
                    constexpr measurement operator/(const measurement& other) const {
                        return measurement(value_ / other.value_, unit_ / other.unit_);
                    }

                    constexpr measurement operator/(const double& val) const {
                        return measurement(value_ / val, unit_);
                    } 

                    constexpr bool operator==(const double& val) const {
                        return (value_ == val) ? true : math::tools::compare_round_equals(value_, val);
                    }

                    constexpr bool operator!=(const double& val) const { return !operator==(val); }

                    constexpr bool operator>(const double& val) const { return value_ > val; }

                    constexpr bool operator<(const double& val) const { return value_ < val; }

                    constexpr bool operator>=(const double& val) const {
                        return (value_ >= val) ? true : operator==(val);
                    }

                    constexpr bool operator<=(const double& val) const {
                        return value_ <= val ? true : operator==(val);
                    }

                    inline bool operator==(const measurement& other) const {
                        return value_equality_check((unit_ == other.units()) ? other.value() : other.value_as(unit_));
                    }

                    inline bool operator!=(const measurement& other) const {
                        return !value_equality_check((unit_ == other.units()) ? other.value() : other.value_as(unit_));
                    }

                    inline bool operator>(const measurement& other) const {
                        return value_ > other.value_as(unit_);
                    }

                    inline bool operator<(const measurement& other) const {
                        return value_ < other.value_as(unit_);
                    }

                    inline bool operator>=(const measurement& other) const {
                        return (value_ > other.value_as(unit_)) ? true : value_equality_check(other.value_as(unit_));
                    }

                    inline bool operator<=(const measurement& other) const {
                        return (value_ < other.value_as(unit_)) ? true : value_equality_check(other.value_as(unit_));
                    }       


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return measurement(math::tools::convert(value_, unit_, newUnits), newUnits);
                    }

                    // convert a unit into its base units
                    constexpr measurement convert_to_base() const {
                        return measurement(value_ * unit_.prefix_.multiplier_, unit_.as_unit());
                    }

                    
                    // =============================================                                                                                         
                    // get methods
                    // =============================================  

                    // get the numerical component of the measurement
                    constexpr double value() const { return value_; }                        
                    
                    // set the numerical component of the measurement
                    constexpr void value(const double& value) { value_ = value; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units);
                    }

                    // get the unit component of a measurement
                    constexpr units::unit units() const { return unit_; }

                    // convert the measurement to a single unit
                    constexpr units::unit as_unit() const { return {value_, unit_ }; }

                    // print the measurement
                    void print() const { 
                        if (unit_.multiplier() != 1) std::cout << std::setprecision((int)-log10(unit_.multiplier())) << value_ << " "; 
                        else std::cout << value_ << " ";
                        unit_.print(); 
                        std::cout << "\n"; 
                    }


                private:

                    // does a numerical equality check on the value accounting for tolerances
                    bool value_equality_check(const double& otherval) const {
                        return (value_ == otherval) ? true : math::tools::compare_round_equals(value_, otherval);
                    } 


            }; // class measurement


            // measurement class with a value and a fixed unit
            class fixed_measurement {

                protected:

                    // =============================================                                                                                         
                    // class members
                    // =============================================  

                    double value_{};  

                    const units::unit unit_;  


                public:

                    // =============================================                                                                                         
                    // constructors & destructor
                    // =============================================  

                    // default constructor
                    fixed_measurement() noexcept {};

                    // constructor from a value and a unit 
                    explicit constexpr fixed_measurement(double val, const units::unit& unit) noexcept :
                        value_(val), unit_(unit) {}

                    // constructor from a measurement
                    constexpr fixed_measurement(const measurement& other) :
                        value_(other.value()), unit_(other.units()) {}

                    // constructor from a fixed measurement
                    fixed_measurement(const fixed_measurement& other) :
                        value_(other.value_), unit_(other.unit_) {}

                    // destructor
                    ~fixed_measurement() = default;


                    // =============================================                                                                                         
                    // operators
                    // =============================================  

                    fixed_measurement& operator=(const measurement& other) noexcept {
                        value_ = (unit_ == other.units()) ? other.value() : other.value_as(unit_);
                        return *this;
                    }

                    fixed_measurement& operator=(const fixed_measurement& other) noexcept {
                        value_ = (unit_ == other.units()) ? other.value() : other.value_as(unit_);
                        return *this;
                    }

                    fixed_measurement& operator=(const double& val) noexcept {
                        value_ = val;
                        return *this;
                    }     

                    fixed_measurement& operator+=(const double& val) {
                        value_ += val;
                        return *this;
                    }

                    fixed_measurement& operator-=(const double& val) {
                        value_ -= val;
                        return *this;
                    }

                    fixed_measurement& operator*=(const double& val) {
                        value_ *= val;
                        return *this;
                    }

                    fixed_measurement& operator/=(const double& val) {
                        value_ /= val;
                        return *this;
                    }

                    constexpr fixed_measurement operator*(const double& val) const {
                        return fixed_measurement(value_ * val, unit_);
                    }

                    constexpr fixed_measurement operator/(const double& val) const {
                        return fixed_measurement(value_ / val, unit_);
                    }

                    inline fixed_measurement operator+(const measurement& other) const {
                        return fixed_measurement(value_ + other.value_as(unit_), unit_);
                    }

                    constexpr fixed_measurement operator+(const double& val) const {
                        return fixed_measurement(value_ + val, unit_);
                    }

                    inline fixed_measurement operator-(const measurement& other) const {
                        return fixed_measurement(value_ - other.value_as(unit_), unit_);
                    }

                    constexpr fixed_measurement operator-(const double& val) const {
                        return fixed_measurement(value_ - val, unit_);
                    }

                    constexpr bool operator==(const double& val) const {
                        return (value_ == val) ? true : math::tools::compare_round_equals(value_, val);
                    }

                    constexpr bool operator!=(const double& val) const { return !operator==(val); }

                    constexpr bool operator>(const double& val) const { return value_ > val; }

                    constexpr bool operator<(const double& val) const { return value_ < val; }

                    constexpr bool operator>=(const double& val) const {
                        return (value_ >= val) ? true : operator==(val);
                    }

                    constexpr bool operator<=(const double& val) const {
                        return value_ <= val ? true : operator==(val);
                    }

                    inline bool operator==(const fixed_measurement& val) const {
                        return operator==((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator!=(const fixed_measurement& val) const {
                        return operator!=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator==(const measurement& val) const {
                        return operator==((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator!=(const measurement& val) const {
                        return operator!=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator>(const measurement& val) const {
                        return operator>((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator<(const measurement& val) const {
                        return operator<((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator>=(const measurement& val) const {
                        return operator>=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }

                    inline bool operator<=(const measurement& val) const {
                        return operator<=((unit_ == val.units()) ? val.value() : val.value_as(unit_));
                    }


                    // =============================================                                                                                         
                    // convert methods
                    // =============================================  

                    // direct conversion operator
                    operator measurement() { return measurement(value_, unit_); }

                    // convert a unit to have a new base
                    measurement convert_to(const units::unit& newUnits) const {
                        return measurement(math::tools::convert(value_, unit_, newUnits), newUnits);
                    }

                    // convert a unit into its base units
                    measurement convert_to_base() const {
                        return measurement(value_ * unit_.multiplier(), units::unit(unit_.base_units()));
                    }


                    // =============================================                                                                                         
                    // get methods
                    // =============================================  

                    // get the numerical component of the measurement
                    constexpr double value() const { return value_; }                        
                    
                    // set the numerical component of the measurement
                    constexpr void value(const double& value) { value_ = value; }

                    // get the numerical value as a particular unit type
                    double value_as(const units::unit& desired_units) const {
                        return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units);
                    }

                    // get the unit component of a measurement
                    constexpr units::unit units() const { return unit_; }

                    // convert the measurement to a single unit
                    constexpr units::unit as_unit() const { return {value_, unit_ }; }

                    // print the measurement
                    void print() const { 
                        if (unit_.multiplier() != 1) std::cout << std::setprecision((int)-log10(unit_.multiplier())) << value_ << " "; 
                        else std::cout << value_ << " ";
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

                    units::unit unit_;       
                

                public: 

                    // =============================================                                                                                         
                    // constructors & destructor 
                    // =============================================  

                    // default constructor
                    uncertain_measurement() noexcept {};

                    // constructor from a value, uncertainty, and unit
                    explicit constexpr uncertain_measurement(const double& val, const double& uncertainty_val, const units::unit& unit) noexcept :
                        value_(val), uncertainty_(uncertainty_val), unit_(unit) {}

                    // constructpr from a value and an unit, assuming the uncertainty is 0
                    explicit constexpr uncertain_measurement(const double& val, const units::unit& unit) noexcept :
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
                    uncertain_measurement simple_product(const uncertain_measurement& other) const {
                        double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                        double nval = value_ * other.value_;
                        return uncertain_measurement(nval, nval * ntol, unit_ * other.units());
                    }

                    // multiply with another measurement equivalent to uncertain_measurement multiplication with 0 uncertainty
                    uncertain_measurement operator*(const measurement& other) const {
                        return uncertain_measurement(value() * other.value(), other.value() * uncertainty(), unit_ * other.units());
                    }

                    uncertain_measurement operator*(const units::unit& other) const {
                        return uncertain_measurement(value_, uncertainty_, unit_ * other);
                    }

                    uncertain_measurement operator*(const double& val) const {
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
                    uncertain_measurement simple_divide(const uncertain_measurement& other) const {
                        double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                        double nval = value_ / other.value_;
                        return uncertain_measurement(nval, nval * ntol, unit_ / other.units());
                    }

                    uncertain_measurement operator/(const measurement& other) const {
                        return uncertain_measurement(static_cast<double>(value() / other.value()),
                            static_cast<double>(uncertainty() / other.value()), unit_ / other.units());
                    }

                    uncertain_measurement operator/(const units::unit& other) const {
                        return uncertain_measurement(value_, uncertainty_, unit_ / other);
                    }

                    uncertain_measurement operator/(const double& val) const {
                        return uncertain_measurement(value_ / val, uncertainty_ / val, unit_);
                    }

                    // compute a unit addition and calculate the new uncertainties using the root sum of squares(rss) method
                    uncertain_measurement operator+(const uncertain_measurement& other) const {
                        auto cval = static_cast<double>(math::tools::convert(other.unit_, unit_));
                        double ntol = math::algebra::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                        return uncertain_measurement(value_ + cval * other.value_, ntol, unit_);
                    }

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

                    uncertain_measurement operator+(const measurement& other) const {
                        auto cval = other.value_as(unit_);
                        return uncertain_measurement(value_ + cval, uncertainty_, unit_);
                    }

                    uncertain_measurement operator-(const measurement& other) const {
                        auto cval = other.value_as(unit_);
                        return uncertain_measurement(value_ - cval, uncertainty_, unit_);
                    }

                    // comparison operators 
                    bool operator==(const measurement& other) const {
                        auto val = other.value_as(unit_);
                        if (uncertainty_ == 0.0F) { return (value_ == val) ? true : math::tools::compare_round_equals(value_, val); }
                        return (val >= (value_ - uncertainty_) && val <= (value_ + uncertainty_));
                    }

                    inline bool operator>(const measurement& other) const {
                        return value_ > other.value_as(unit_);
                    }

                    inline bool operator<(const measurement& other) const {
                        return value_ < other.value_as(unit_);
                    }

                    bool operator>=(const measurement& other) const {
                        auto val = other.value_as(unit_);
                        return (value() >= val) ? true : operator==(measurement(val, unit_));
                    }

                    bool operator<=(const measurement& other) const {
                        auto val = other.value_as(unit_);
                        return (value() <= val) ? true : operator==(measurement(val, unit_));
                    }

                    inline bool operator!=(const measurement& other) const {
                        return !operator==(other);
                    }

                    bool operator==(const uncertain_measurement& other) const {
                        auto zval = simple_subtract(other);
                        return (zval == measurement(0.0, unit_));
                    }

                    inline bool operator>(const uncertain_measurement& other) const {
                        return value_ > other.value_as(unit_);
                    }

                    inline bool operator<(const uncertain_measurement& other) const {
                        return value_ < other.value_as(unit_);
                    }

                    bool operator>=(const uncertain_measurement& other) const {
                        auto zval = simple_subtract(other);
                        return (zval.value_ >= 0.0F) ? true :
                                                    (zval == measurement(0.0, unit_));
                    }

                    bool operator<=(const uncertain_measurement& other) const {
                        auto zval = simple_subtract(other);
                        return (zval.value_ <= 0.0F) ? true : (zval == measurement(0.0, unit_));
                    }

                    inline bool operator!=(const uncertain_measurement& other) const {
                        return !operator==(other);
                    }


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
                    constexpr units::unit units() const { return unit_; }

                    // get the numerical value as a particular unit type
                    inline double value_as(const units::unit& desired_units) const {
                        return (unit_ == desired_units) ? value_ : math::tools::convert(value_, unit_, desired_units);
                    }

                    // get the numerical value of the uncertainty as a particular unit
                    inline double uncertainty_as(const units::unit& desired_units) const {
                        return (unit_ == desired_units) ? uncertainty_ : math::tools::convert(uncertainty_, unit_, desired_units);
                    }

                    // convert a unit to have a new base
                    uncertain_measurement convert_to(const units::unit& newUnits) const {
                        auto cval = math::tools::convert(unit_, newUnits);
                        return uncertain_measurement(cval * value_, uncertainty_ * cval, newUnits);
                    }

                    // print the uncertain measurement
                    void print() const { 
                        if (unit_.multiplier() != 1) std::cout << std::setprecision((int)-log10(unit_.multiplier())) << value() << " ± " << uncertainty() << " "; 
                        else std::cout << value() << " ± " << uncertainty() << " "; 
                        unit_.print(); 
                        std::cout << "\n";
                    }


            }; // class uncertain_measurement


        } // namespace measurements


    } // namespace physics


    namespace math {

        // namespace defining some usefull operation
        namespace algebra {

            // generate a unit which is an integer power of another
            constexpr physics::units::unit pow(const physics::units::unit& u, const int& power) { return u.pow(power); }

            // generate the square of an physics::units::unit
            constexpr physics::units::unit square(const physics::units::unit& u) { return u.pow(2); }

            // generate the cube of an physics::units::unit
            constexpr physics::units::unit cube(const physics::units::unit& u) { return u.pow(3); }

            // generate the root of an physics::units::unit
            physics::units::unit root(const physics::units::unit& un, const int& power) {
                if (power == 0) { return physics::units::special::one; }
                if (un.multiplier() < 0.0 && power % 2 == 0) { return physics::units::special::invalid; }
                return physics::units::unit(un.base_units().root(power), root(un.multiplier(), power));
            }

            // generate the square root of an physics::units::unit
            inline physics::units::unit sqrt(const physics::units::unit& u) { return sqrt(u); }

            // generate the cubic root of an physics::units::unit
            inline physics::units::unit cbrt(const physics::units::unit& u) { return cbrt(u); }

            constexpr physics::measurements::measurement operator*(const double& val, const physics::measurements::measurement& meas) { return meas * val; }

            constexpr physics::measurements::measurement operator*(const double& val, const physics::units::unit& unit_base) { return physics::measurements::measurement(val, unit_base); }

            constexpr physics::measurements::measurement operator*(const physics::units::unit& unit_base, const double& val) { return physics::measurements::measurement(val, unit_base); }

            constexpr physics::measurements::fixed_measurement operator*(const double& v1, const physics::measurements::fixed_measurement& v2) { return physics::measurements::fixed_measurement(v1 * v2.value(), v2.units()); }

            constexpr physics::measurements::fixed_measurement operator*(const physics::measurements::fixed_measurement& v2, const double& v1) { return physics::measurements::fixed_measurement(v1 * v2.value(), v2.units()); }

            inline physics::measurements::uncertain_measurement operator*(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) { return v2.operator*(v1); }

            inline physics::measurements::uncertain_measurement operator*(const double& v1, const physics::measurements::uncertain_measurement& v2) { return v2.operator*(v1); }

            inline physics::measurements::uncertain_measurement operator*(const physics::measurements::uncertain_measurement& v2, const double& v1) { return v2.operator*(v1); }


            constexpr physics::measurements::measurement operator/(const double& val, const physics::measurements::measurement& meas) {
                return physics::measurements::measurement(val / meas.value(), meas.units().inv());
            }

            constexpr physics::measurements::measurement operator/(const physics::measurements::measurement& meas, const double& val) {
                return physics::measurements::measurement(val / meas.value(), meas.units().inv());
            }
            
            constexpr physics::measurements::measurement operator/(const double& val, const physics::units::unit& unit_base) {
                return physics::measurements::measurement(val, unit_base.inv());
            }
            
            constexpr physics::measurements::measurement operator/(const physics::units::unit& unit_base, const double& val) {
                return physics::measurements::measurement(1.0 / val, unit_base);
            }

            constexpr physics::measurements::fixed_measurement operator/(const double& v1, const physics::measurements::fixed_measurement& v2) {
                return physics::measurements::fixed_measurement(v1 / v2.value(), v2.units().inv());
            }                
            
            constexpr physics::measurements::fixed_measurement operator/(const physics::measurements::fixed_measurement& v2, const double& v1) {
                return physics::measurements::fixed_measurement(v1 / v2.value(), v2.units().inv());
            }

            physics::measurements::uncertain_measurement operator/(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
                double ntol = v2.uncertainty() / v2.value();
                double nval = v1.value() / v2.value();
                return physics::measurements::uncertain_measurement(nval, nval * ntol, v1.units() / v2.units());
            }

            physics::measurements::uncertain_measurement operator/(const double& v1, const physics::measurements::uncertain_measurement& v2) {
                double ntol = v2.uncertainty() / v2.value();
                double nval = v1 / v2.value();
                return physics::measurements::uncertain_measurement(nval, nval * ntol, v2.units().inv());
            }

            physics::measurements::uncertain_measurement operator/(const int& v1, const physics::measurements::uncertain_measurement& v2) {
                double ntol = v2.uncertainty() / v2.value();
                double nval = static_cast<double>(v1) / v2.value();
                return physics::measurements::uncertain_measurement(nval, nval * ntol, v2.units().inv());
            }


            constexpr physics::measurements::fixed_measurement operator+(const double& v1, const physics::measurements::fixed_measurement& v2) {
                return physics::measurements::fixed_measurement(v1 + v2.value(), v2.units());
            }

            constexpr physics::measurements::fixed_measurement operator-(const double& v1, const physics::measurements::fixed_measurement& v2) {
                return physics::measurements::fixed_measurement(v1 - v2.value(), v2.units());
            }
            
            physics::measurements::uncertain_measurement operator+(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
                double cval = math::tools::convert(v2.units(), v1.units());
                double ntol = v2.uncertainty() * cval;
                return physics::measurements::uncertain_measurement(v1.value() + cval * v2.value(), ntol, v1.units());
            }

            physics::measurements::uncertain_measurement operator-(const physics::measurements::measurement& v1, const physics::measurements::uncertain_measurement& v2) {
                double cval = math::tools::convert(v2.units(), v1.units());
                double ntol = v2.uncertainty() * cval;
                return physics::measurements::uncertain_measurement(v1.value() - cval * v2.value(), ntol, v1.units());
            }


            constexpr physics::measurements::measurement pow(const physics::measurements::measurement& meas, const int& power) {
                return physics::measurements::measurement{ std::pow(meas.value(), power), meas.units().pow(power) };
            }

            constexpr physics::measurements::fixed_measurement pow(const physics::measurements::fixed_measurement& meas, const int& power) {
                return physics::measurements::fixed_measurement{std::pow(meas.value(), power), meas.units().pow(power)};
            }

            physics::measurements::uncertain_measurement pow(const physics::measurements::uncertain_measurement& meas, const int& power) {
                auto new_value = std::pow(meas.value(), power);
                auto new_tol = ((power >= 0) ? power : -power) * new_value * meas.uncertainty() / meas.value();
                return physics::measurements::uncertain_measurement(new_value, new_tol, meas.units().pow(power));                    
            }

            constexpr physics::measurements::measurement square(const physics::measurements::measurement& meas) { return square(meas); }
            
            constexpr physics::measurements::fixed_measurement square(const physics::measurements::fixed_measurement& meas) { return square(meas); }

            constexpr physics::measurements::uncertain_measurement square(const physics::measurements::uncertain_measurement& meas) { return square(meas); }

            constexpr physics::measurements::measurement cube(const physics::measurements::measurement& meas) { return cube(meas); }
            
            constexpr physics::measurements::fixed_measurement cube(const physics::measurements::fixed_measurement& meas) { return cube(meas); }
            
            constexpr physics::measurements::uncertain_measurement cube(const physics::measurements::uncertain_measurement& meas) { return cube(meas); }


            physics::measurements::measurement root(const physics::measurements::measurement& meas, const int& power) {
                return physics::measurements::measurement(root(meas.value(), power), root(meas.units(), power));
            }

            physics::measurements::fixed_measurement root(const physics::measurements::fixed_measurement& meas, const int& power) {
                return physics::measurements::fixed_measurement(root(meas.value(), power), root(meas.units(), power));
            }

            physics::measurements::uncertain_measurement root(const physics::measurements::uncertain_measurement& um, const int& power) {
                auto new_value = root(um.value(), power);
                auto new_tol = new_value * um.uncertainty() / (static_cast<double>((power >= 0) ? power : -power) * um.value());
                return physics::measurements::uncertain_measurement(new_value, new_tol, root(um.units(), power));
            }

            inline physics::measurements::measurement sqrt(const physics::measurements::measurement& meas) { return sqrt(meas); }
            
            inline physics::measurements::fixed_measurement sqrt(const physics::measurements::fixed_measurement& meas) { return sqrt(meas); }

            inline physics::measurements::uncertain_measurement sqrt(const physics::measurements::uncertain_measurement& meas) { return sqrt(meas); }            

            inline physics::measurements::measurement cbrt(const physics::measurements::measurement& meas) { return cbrt(meas); }

            inline physics::measurements::fixed_measurement cbrt(const physics::measurements::fixed_measurement& meas) { return cbrt(meas); }

            inline physics::measurements::uncertain_measurement cbrt(const physics::measurements::uncertain_measurement& meas) { return cbrt(meas); }


            constexpr bool operator==(const double& val, const physics::measurements::fixed_measurement& v2) {
                return v2 == val;
            }

            constexpr bool operator!=(const double& val, const physics::measurements::fixed_measurement& v2) {
                return v2 != val;
            }

            constexpr bool operator>(const double& val, const physics::measurements::fixed_measurement& v2) {
                return val > v2.value();
            }

            constexpr bool operator<(const double& val, const physics::measurements::fixed_measurement& v2) {
                return val < v2.value();
            }

            constexpr bool operator>=(const double& val, const physics::measurements::fixed_measurement& v2) {
                return (val >= v2.value()) ? true : (v2 == val);
            }

            constexpr bool operator<=(const double& val, const physics::measurements::fixed_measurement& v2) {
                return (val <= v2.value()) ? true : (v2 == val);
            }

            inline bool operator==(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return v2 == other;
            }

            inline bool operator!=(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return v2 != other;
            }

            constexpr bool operator>(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return other.value() > v2.value();
            }

            constexpr bool operator<(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return other.value() < v2.value();
            }

            constexpr bool operator>=(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return (other > v2) ? true : (v2 == other);
            }

            constexpr bool operator<=(const physics::measurements::measurement& other, const physics::measurements::uncertain_measurement& v2) {
                return (other < v2) ? true : (v2 == other);
            }


        } // namespace algebra


    } // namespace math


} // namespace physim
