#include "../include/spring.h"

#include <cmath>

atg_scs::Spring::Spring() {
    m_restLength = 1.0;
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
    if (m_body1 == nullptr || m_body2 == nullptr) return;

    double x1, y1;
    double x2, y2;

    double v_x1 = 0, v_y1 = 0;
    double v_x2 = 0, v_y2 = 0;

    if (m_body1->index != -1) {
        state->localToWorld(m_p1_x, m_p1_y, &x1, &y1, m_body1->index);
        state->velocityAtPoint(m_p1_x, m_p1_y, &v_x1, &v_y1, m_body1->index);
    }
    else {
        m_body1->localToWorld(m_p1_x, m_p1_y, &x1, &y1);
    }

    if (m_body2->index != -1) {
        state->localToWorld(m_p2_x, m_p2_y, &x2, &y2, m_body2->index);
        state->velocityAtPoint(m_p2_x, m_p2_y, &v_x2, &v_y2, m_body2->index);
    }
    else {
        m_body2->localToWorld(m_p2_x, m_p2_y, &x2, &y2);
    }

    double dx = x2 - x1;
    double dy = y2 - y1;

    const double l = std::sqrt(dx * dx + dy * dy);

    if (std::abs(l) >= 1E-2) {
        dx /= l;
        dy /= l;
    }
    else {
        dx = 0.0;
        dy = 0.0;
    }

    const double rel_v_x = (v_x2 - v_x1);
    const double rel_v_y = (v_y2 - v_y1);

    const double v = dx * rel_v_x + dy * rel_v_y;
    const double x = l - m_restLength;

    state->applyForce(
        m_p1_x,
        m_p1_y,
        dx * x * m_ks + rel_v_x * m_kd,
        dy * x * m_ks + rel_v_y * m_kd,
        m_body1->index
    );

    state->applyForce(
        m_p2_x,
        m_p2_y,
        -dx * x * m_ks - rel_v_x * m_kd,
        -dy * x * m_ks - rel_v_y * m_kd,
        m_body2->index
    );
}

void atg_scs::Spring::getEnds(double *x_1, double *y_1, double *x_2, double *y_2) {
    if (m_body1 == nullptr || m_body2 == nullptr) return;

    m_body1->localToWorld(m_p1_x, m_p1_y, x_1, y_1);
    m_body2->localToWorld(m_p2_x, m_p2_y, x_2, y_2);
}

double atg_scs::Spring::energy() const {
    if (m_body1 == nullptr || m_body2 == nullptr) return 0;

    double x1, y1;
    double x2, y2;

    m_body1->localToWorld(m_p1_x, m_p1_y, &x1, &y1);
    m_body2->localToWorld(m_p2_x, m_p2_y, &x2, &y2);

    const double dx = x2 - x1;
    const double dy = y2 - y1;

    const double l = std::sqrt(dx * dx + dy * dy);

    return 0.5 * m_ks * (l - m_restLength) * (l - m_restLength);
}
