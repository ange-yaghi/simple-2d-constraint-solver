#include "../include/spring.h"

#include <cmath>

atg_scs::Spring::Spring() {
    m_restLength = 0;
    m_ks = 0;
    m_kd = 0;

    m_p1_x = m_p1_y = 0;
    m_p2_x = m_p2_y = 0;

    m_body1 = m_body2 = nullptr;
}

atg_scs::Spring::~Spring() {
    /* void */
}

void atg_scs::Spring::apply(SystemState *state) {
    double x1, y1;
    double x2, y2;

    double v_x1, v_y1;
    double v_x2, v_y2;

    state->localToWorld(m_p1_x, m_p1_y, &x1, &y1, m_body1->index);
    state->localToWorld(m_p2_x, m_p2_y, &x2, &y2, m_body2->index);

    state->velocityAtPoint(m_p1_x, m_p1_y, &v_x1, &v_y1, m_body1->index);
    state->velocityAtPoint(m_p2_x, m_p2_y, &v_x2, &v_y2, m_body2->index);
    
    double dx = x2 - x1;
    double dy = y2 - y1;

    const double l = std::sqrt(dx * dx + dy * dy);

    dx /= l;
    dy /= l;

    const double rel_v_x = (v_x2 - v_x1);
    const double rel_v_y = (v_y2 - v_y1);

    const double v = dx * rel_v_x + dy * rel_v_y;

    state->applyForce(
        m_p1_x,
        m_p1_y,
        dx * ((l - m_restLength) * m_ks + v * m_kd),
        dy * ((l - m_restLength) * m_ks + v * m_kd),
        m_body1->index
    );

    state->applyForce(
        m_p2_x,
        m_p2_y,
        -dx * ((l - m_restLength) * m_ks + v * m_kd),
        -dy * ((l - m_restLength) * m_ks + v * m_kd),
        m_body2->index
    );
}


void atg_scs::Spring::getEnds(double *x_1, double *y_1, double *x_2, double *y_2) {
    m_body1->localToWorld(m_p1_x, m_p1_y, x_1, y_1);
    m_body2->localToWorld(m_p2_x, m_p2_y, x_2, y_2);
}
