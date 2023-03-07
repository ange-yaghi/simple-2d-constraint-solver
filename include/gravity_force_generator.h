#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_GRAVITY_FORCE_GENERATOR_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_GRAVITY_FORCE_GENERATOR_H

#define SUN_GRAVITY (274)
#define EARTH_GRAVITY (9.81)
//reference: https://nssdc.gsfc.nasa.gov/planetary/factsheet/planet_table_ratio.html
#define MERCURY_GRAVITY (0.378*EARTH_GRAVITY)
#define VENUS_GRAVITY (0.907*EARTH_GRAVITY)
#define MOON_GRAVITY (0.166*EARTH_GRAVITY)
#define MARS_GRAVITY (0.377*EARTH_GRAVITY)
#define JUPITER_GRAVITY (2.36*EARTH_GRAVITY)
#define SATURN_GRAVITY (0.916*EARTH_GRAVITY)
#define URANUS_GRAVITY (0.889*EARTH_GRAVITY)
#define NEPTUNE_GRAVITY (1.12*EARTH_GRAVITY)
#define PLUTO_GRAVITY (0.071*EARTH_GRAVITY)

#include "force_generator.h"

#include "rigid_body.h"

namespace atg_scs {
    class GravityForceGenerator : public ForceGenerator {
        public:
            GravityForceGenerator();
            virtual ~GravityForceGenerator();

            virtual void apply(SystemState *state);

            double m_g;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_GRAVITY_FORCE_GENERATOR_H */
