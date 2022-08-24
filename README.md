# physim
c++ header only namespace for computational physics

namespace physim {

    namespace physics {

        namespace constants {}

        namespace units {
            class unit_data {}; 
            class unit_prefix {}; 
            class unit : public unit_data, public unit_prefix {}; 
            namespace defined {}
        }

        namespace measurements {
            class measurement {};
            class fixed_measurement {};
            class uncertain_measurement {}; 
        }

        namespace position {
            class coordinate {}; 
            class position {}; 
            class velocity {}; 
        }

        namespace objects {
            class mass {}; 
            class charge {}; 
        }

        namespace time {}

    }

    namespace math {

        namespace functions {}
        namespace statistics {}
        class tools::random_generator {};
        class integral {}; 

    }