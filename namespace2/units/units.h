
// author:          Lorenzo Liuzzo
// email:           lorenzoliuzzo@outlook.com
// description:     physics(namespace) containing the basic tools for computational physics. 
// last updated:    27/07/2022

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip> 
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>


namespace physim {

    namespace physics {

        namespace tools {

            // namespace defining some usefull constants
            namespace constants {

                constexpr int32_t max_neg(uint32_t n_bits) { return -(int32_t(1U << (n_bits - 1))); }
                constexpr double infinity = std::numeric_limits<double>::infinity();
                constexpr double invalid_conversion = std::numeric_limits<double>::signaling_NaN();

                constexpr double pi = 3.14159265358979323846;
                constexpr double e = 2.7182818284590452353603;

            } // namespace constants


            // namespace defining some usefull operation
            namespace op {

                // generate the square power of a value
                template<typename T> constexpr T square(const T& a) { return a * a; }

                // generate the cubic power of a value
                template<typename T> constexpr T cube(const T& a) { return a * a * a; }

                // generate small integer powers of a value(1, 0, -1)
                template<typename T> constexpr T pow_small(const T& val, const int& power) { 
                    return (power == 1) ? val : ((power == -1) ? T(1.0) / val : T(1.0));
                }

                // generate an integer power of a number
                template<typename T> constexpr T pow(T val, const int& power) {
                    return (power > 1) ? square(pow(val, power / 2)) * (power % 2 == 0 ? T(1.0) : val) :
                        (power < -1) ? T(1.0) / (square(pow(val, (-power) / 2)) * ((-power) % 2 == 0 ? T(1.0) : val)) :
                        pow_small(val, power);
                }

                // generate root power of a value
                template<typename X> X root(X value, const int& power) {
                    switch (power) {
                        case 0: return X{1.0};
                        case 1: return value;
                        case -1: return X{1.0} / value;
                        case 2: 
                            if (value < X{0.0}) { return constants::invalid_conversion; }
                            else return std::sqrt(value);
                        case -2:
                            if (value < X{0.0}) { return constants::invalid_conversion; }
                            else return std::sqrt(X{1.0} / value);
                        case 3: return std::cbrt(value);
                        case -3: return std::cbrt(X{1.0} / value);
                        case 4: 
                            if (value < X{0.0}) { return constants::invalid_conversion; }
                            else return std::sqrt(std::sqrt(value));
                        case -4: 
                            if (value < X{0.0}) { return constants::invalid_conversion; }
                            else return std::sqrt(std::sqrt(X{1.0} / value));
                        default:
                            if (value < X{0.0} && power % 2 == 0) { return constants::invalid_conversion; }
                            else return std::pow(value, X{1.0} / static_cast<X>(power));
                    }
                }
                
                // generate the square root power of a value
                template <typename X> X sqrt(X x) { return root(x, 2); }

                // generate the cubic root power of a value
                template <typename X> X cbrt(X x) { return root(x, 3); }

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

            } // namespace op

            
            // namespace defining the units of measurements
            namespace units {
                
                // number of bits used for encoding base unit exponents 
                namespace bitwidth {

                    constexpr uint32_t base_size = sizeof(uint32_t) == 8 ? 8 : 4;
                    constexpr uint32_t meter{(base_size == 8) ? 8 : 4};
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

                        signed int meter_ : bitwidth::meter;
                        signed int second_ : bitwidth::second;  
                        signed int kilogram_ : bitwidth::kilogram;
                        signed int ampere_ : bitwidth::ampere;
                        signed int candela_ : bitwidth::candela;  
                        signed int kelvin_ : bitwidth::kelvin;
                        signed int mole_ : bitwidth::mole;


                    public:

                        // the seven SI base units 
                        enum base {
                            Meter = 0,
                            Second = 1,
                            Kilogram = 2,
                            Ampere = 3,
                            Kelvin = 4,
                            Mole = 5,
                            Candela = 6
                        };

                        static constexpr uint32_t bits[7] = { 
                            bitwidth::meter, bitwidth::second, 
                            bitwidth::kilogram, bitwidth::ampere, 
                            bitwidth::kelvin, bitwidth::mole, 
                            bitwidth::candela
                        };  


                        // =============================================
                        // constructors
                        // ============================================= 
                        
                        // constructor from powers
                        constexpr unit_data(const int& meters, const int& seconds, const int& kilograms, 
                                            const int& amperes, const int& kelvins,  const int& moles, const int& candelas) :
                            meter_(meters), second_(seconds), kilogram_(kilograms), ampere_(amperes),
                            kelvin_(kelvins), mole_(moles), candela_(candelas) {};

                        explicit constexpr unit_data(std::nullptr_t) :
                            meter_(constants::max_neg(bitwidth::meter)), second_(constants::max_neg(bitwidth::second)), 
                            kilogram_(constants::max_neg(bitwidth::kilogram)), ampere_(constants::max_neg(bitwidth::ampere)),
                            kelvin_(constants::max_neg(bitwidth::kelvin)), mole_(constants::max_neg(bitwidth::mole)), 
                            candela_(constants::max_neg(bitwidth::candela)) {}


                        // =============================================
                        // operators
                        // ============================================= 

                        // perform a multiply operation by adding the powers together
                        constexpr unit_data operator*(const unit_data& other) const {
                            return { 
                                meter_ + other.meter_,
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
                                meter_ - other.meter_,
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
                                -meter_,
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
                                meter_ * power,
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
                            if (has_valid_root(power)) return unit_data( meter_ / power,
                                                                        second_ / power,
                                                                        kilogram_ / power,
                                                                        ampere_ / power,
                                                                        kelvin_ / power,
                                                                        mole_ / power,
                                                                        candela_ / power);
                            else exit(-11);
                        }
                        
                        
                        // =============================================
                        // check methods
                        // ============================================= 

                        // comparison operators
                        constexpr bool operator==(const unit_data& other) const {
                            return meter_ == other.meter_ && second_ == other.second_ &&
                                kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                                candela_ == other.candela_ && kelvin_ == other.kelvin_ && mole_ == other.mole_;         
                        }

                        constexpr bool operator!=(const unit_data& other) const {
                            return !(*this == other);
                        }

                        // check if the units have the same base unit 
                        constexpr bool has_same_base(const unit_data& other) const {
                            return *this == other;
                        }

                        // check if the unit is empty
                        constexpr bool empty() const {
                            return meter_ == 0 && second_ == 0 && kilogram_ == 0 &&
                                ampere_ == 0 && candela_ == 0 && kelvin_ == 0 && mole_ == 0;
                        }
                        

                        // =============================================
                        // get methods
                        // =============================================
                        
                        constexpr int meter() const { return meter_; }

                        constexpr int second() const { return second_; }

                        constexpr int kg() const { return kilogram_; }
                                            
                        constexpr int ampere() const { return ampere_; }
                        
                        constexpr int kelvin() const { return kelvin_; }
                        
                        constexpr int mole() const { return mole_; }
                        
                        constexpr int candela() const { return candela_; }

                        constexpr int unit_type_count() const {
                            return ((meter_ != 0) ? 1 : 0) + ((second_ != 0) ? 1 : 0) +
                                ((kilogram_ != 0) ? 1 : 0) + ((ampere_ != 0) ? 1 : 0) +
                                ((kelvin_ != 0) ? 1 : 0) + ((mole_ != 0) ? 1 : 0) + 
                                ((candela_ != 0) ? 1 : 0);
                        }

                        constexpr unit_data base_units() const { return *this; }

                        constexpr void print() const {
                            if (meter_ != 0 && meter_ != 1) std::cout << "m^" << meter_ << "\n"; 
                            if (meter_ == 1) std::cout << "m" << "\n";
                            if (second_ != 0 && second_ != 1) std::cout << "s^" << second_ << "\n"; 
                            if (second_ == 1) std::cout << "s" << "\n"; 
                            if (kilogram_ != 0 && kilogram_ != 1) std::cout << "kg^" << kilogram_ << "\n"; 
                            if (kilogram_ == 1) std::cout << "kg" << "\n"; 
                            if (ampere_ != 0 && ampere_ != 1) std::cout << "A^" << ampere_ << "\n"; 
                            if (ampere_ == 1) std::cout << "A" << "\n"; 
                            if (kelvin_ != 0 && kelvin_ != 1) std::cout << "K^" << kelvin_ << "\n"; 
                            if (kelvin_ == 1) std::cout << "K" << "\n";
                            if (mole_ != 0 && mole_ != 1) std::cout << "mol^" << mole_ << "\n"; 
                            if (mole_ == 1) std::cout << "mol"  << "\n"; 
                            if (candela_ != 0 && candela_ != 1) std::cout << "cd^" << candela_ << "\n"; 
                            if (candela_ == 1) std::cout << "cd" << "\n"; 
                        }


                    private: 

                        // check if the base_unit has a valid root
                        constexpr bool has_valid_root(const int& power) const {
                            return meter_ % power == 0 && second_ % power == 0 &&
                                kilogram_ % power == 0 && ampere_ % power == 0 &&
                                candela_ % power == 0 && kelvin_ % power == 0 &&
                                mole_ % power == 0;
                        }      
                        
                        // to handle a few weird operations that operate on square_root Hz
                        constexpr int root_Hertz_modifier(const int& power) const {
                            return (second_ * power == 0 || 
                                    power % 2 != 0) ? 0 : (power / 2) * ((second_ < 0) || 
                                    (power < 0) ? 9 : -9);
                        }    

                }; // class unit_base


                // class defining a basic unit module 
                class unit_prefix {

                    protected: 

                        // =============================================
                        // class members
                        // ============================================= 
                        
                        double multiplier_{1.0};  

                        const char prefix_ = '-'; 


                    public: 

                        // =============================================
                        // constructors
                        // ============================================= 

                        // default constructor
                        constexpr unit_prefix() noexcept {};

                        explicit constexpr unit_prefix(const double& mult) noexcept : 
                            multiplier_{mult} {} 

                        explicit constexpr unit_prefix(const double& mult, const char& pref) noexcept : 
                            multiplier_{mult}, prefix_(pref) {} 

                        explicit constexpr unit_prefix(const char& pref, const double& mult) noexcept : 
                            multiplier_{mult}, prefix_(pref) {} 
                            

                        // =============================================
                        // operators
                        // ============================================= 

                        // perform a multiply operation by adding the powers together
                        constexpr unit_prefix operator*(const unit_prefix& other) const {
                            return unit_prefix(multiplier_ + other.multiplier(), prefix_ + other.prefix_);
                        }

                        // perform a division operation by subtract the powers together
                        constexpr unit_prefix operator/(const unit_prefix& other) const {
                            return unit_prefix(multiplier_ / other.multiplier(), prefix_);
                        }

                        // invert the unit
                        constexpr unit_prefix inv() const {
                            return unit_prefix(1 / multiplier_, prefix_);
                        }

                        // take a unit_prefix to some power
                        constexpr unit_prefix pow(const int& power) const { 
                            return unit_prefix(op::pow(multiplier_, power), prefix_);
                        }
                        
                        // take some root of a unit_prefix
                        unit_prefix root(const int& power) const {
                            return unit_prefix(op::root(multiplier_, power), prefix_); 
                        }
                        
                        
                        // =============================================
                        // check methods
                        // ============================================= 

                        // comparison operators
                        constexpr bool operator==(const unit_prefix& other) const {
                            if (multiplier_ == other.multiplier_) return true;    
                            else return op::compare_round_equals(multiplier_, other.multiplier_);
                        }

                        constexpr bool operator!=(const unit_prefix& other) const {
                            return !(*this == other);
                        }

                        // check if the units have the same prefix unit 
                        constexpr bool has_same_prefix(const unit_prefix& other) const {
                            return *this == other;
                        }

                        // check if the multiplier is nan
                        inline bool is_nan(const unit_prefix& pref) {
                            return std::isnan(pref.multiplier_);
                        }

                        // checks that the multiplier is finite
                        inline bool is_finite(const unit_prefix& pref) {
                            return std::isfinite(pref.multiplier_);
                        }

                        // check if the multiplier is infinite
                        inline bool is_inf(const unit_prefix& pref) {
                            return std::isinf(pref.multiplier_);
                        }


                        // =============================================
                        // get methods
                        // ============================================= 
                        
                        constexpr double multiplier() const { return multiplier_; }

                        constexpr unit_prefix prefix() const { return *this; }

                        constexpr void print() const { 
                            if (prefix_ != '-') std::cout << prefix_; 
                            else {
                                if (multiplier_ == 1.0) std::cout << ""; 
                                if (multiplier_ == 1e-1) std::cout << "d"; 
                                if (multiplier_ == 1e-2) std::cout << "c"; 
                                if (multiplier_ == 1e-3) std::cout << "m"; 
                                if (multiplier_ == 1e-6) std::cout << "u"; 
                                if (multiplier_ == 1e-9) std::cout << "n"; 
                                if (multiplier_ == 1e-12) std::cout << "p"; 
                                if (multiplier_ == 1e-15) std::cout << "f"; 
                                if (multiplier_ == 1e-18) std::cout << "a"; 
                                if (multiplier_ == 1e-21) std::cout << "z"; 
                                if (multiplier_ == 1e-24) std::cout << "y"; 
                                if (multiplier_ == 1e2) std::cout << "h"; 
                                if (multiplier_ == 1e3) std::cout << "k"; 
                                if (multiplier_ == 1e6) std::cout << "M"; 
                                if (multiplier_ == 1e9) std::cout << "G"; 
                                if (multiplier_ == 1e12) std::cout << "T"; 
                                if (multiplier_ == 1e15) std::cout << "P"; 
                                if (multiplier_ == 1e18) std::cout << "E"; 
                                if (multiplier_ == 1e21) std::cout << "Z"; 
                                if (multiplier_ == 1e24) std::cout << "Y"; 
                            }
                        }

                }; // class unit_prefix

            } // namespace units

        } // namespace tools

    } // namespace physics

} // namespace physim