#include "../include/gravity_force_generator.h"

//may need improvement - reference: https://nssdc.gsfc.nasa.gov/planetary/factsheet/planet_table_ratio.html
#define SUN_GRAVITY (274)
#define EARTH_GRAVITY (9.81)
#define MERCURY_GRAVITY (0.378*EARTH_GRAVITY)
#define VENUS_GRAVITY (0.907*EARTH_GRAVITY)
#define MOON_GRAVITY (0.166*EARTH_GRAVITY)
#define MARS_GRAVITY (0.377*EARTH_GRAVITY)
#define JUPITER_GRAVITY (2.36*EARTH_GRAVITY)
#define SATURN_GRAVITY (0.916*EARTH_GRAVITY)
#define URANUS_GRAVITY (0.889*EARTH_GRAVITY)
#define NEPTUNE_GRAVITY (1.12*EARTH_GRAVITY)
#define PLUTO_GRAVITY (0.071*EARTH_GRAVITY)

atg_scs::GravityForceGenerator::GravityForceGenerator() {
    m_g = EARTH_GRAVITY;
}

atg_scs::GravityForceGenerator::~GravityForceGenerator() {
    /* void */
}

void atg_scs::GravityForceGenerator::apply(SystemState *state) {
    const int n = state->n;

    for (int i = 0; i < n; ++i) {
        state->f_y[i] += -state->m[i] * m_g;
    }
}
