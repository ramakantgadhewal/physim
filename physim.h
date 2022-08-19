
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physim(namespace) containing the basic tools for computational physics. 
// last updated:    19/08/2022

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
                std::cerr << "Error! It was not possible to open the selected file" << std::endl; 
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
                
                void print() { std::cout << "elapsed time = " << elapsed_time() << std::endl; }

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

            // generate an integer power of a value
            template<typename T> 
            constexpr T pow(const T& value, const int& power) {
                return (power > 1) ? square(pow(value, power / 2)) * (power % 2 == 0 ? T(1.0) : value) :
                    (power < -1) ? T(1.0) / (square(pow(value, (-power) / 2)) * ((-power) % 2 == 0 ? T(1.0) : value)) :
                    pow_small(value, power);
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
                for (auto x : v) accu += pow(x - expected_value, 2) / sd(v); 
                return accu; 
            }           

            // chi squared reduced of an std::vector
            template <typename K>
            double chi_sq_r(const std::vector<K>& v, const double& expected_value, const unsigned int& gdl) {
                double accu{}; 
                for (auto x : v) accu += pow(x - expected_value, 2) / sd(v); 
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
                        std::cout << "f (" << x << ") = " << std::fixed << std::setprecision((int)-log10(precision)) << eval(x) << std::endl; 
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
                        : m_a{a}, m_b{b}, m_c{c}, m_delta{pow(m_b, 2) - 4 * m_a * m_c} {}
                    
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
                        else std::cout << std::endl << "There are not real solutions..." << std::endl << std::endl; // note: create a complex class then come back here
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
                        if (m_a != 1) std::cout << m_a; 
                        std::cout << "x^2 "; 
                        if (m_b != 0) {
                            if (m_b != 1) std::cout << "+ " << m_b << "x "; 
                            else std::cout << "+ x ";                       
                        }
                        if (m_c != 0) { std::cout << "+ " << m_c << std::endl; }
                        else std::cout << std::endl;
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
                        if (m_d != 0) std::cout << " + " << m_d << std::endl;
                        else std::cout << std::endl;  
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
                        std::cout << "x^(1/2)" << std::endl;
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
                        std::cout << "x^(1/3)" << std::endl;
                    }        
                

            }; // class cubic_root


            // exponential function
            class exponential : public function_base {
            
                private:
            
                    // =============================================
                    // class members
                    // =============================================
                    
                    // y = c1 * (base ^ (c2 * x))
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
            
                    constexpr double eval(const double& x) const override { return m_c1 * math::algebra::pow(m_base, m_c2 * x); }
                            
                
                    // =============================================
                    // print methods
                    // =============================================
            
                    void print_equation() const override {
                        std::cout << " y = "; 
                        if (m_c1 != 1) std::cout << m_c1;
                        if (m_base != math::constants::e) std::cout << m_base << "^"; 
                        else std::cout << "e^(";
                        if (m_c2 != 1) std::cout << m_c2; 
                        std::cout << "x)" << std::endl; 
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
                        if (m_c2 != 1) std::cout << "(" << m_c2 << "x)" << std::endl; 
                        else std::cout << "(x)" << std::endl;
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
                        std::cout << "x)" << std::endl;
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
                        std::cout << "x)" << std::endl;
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
                        std::cout << "x)" << std::endl;
                    }


            }; // class tangent


        } // namespace functions
        

    } // namespace math

} // namespace physim
