#include "lift_calculator.h"
#include <cmath>

LiftCalculator::LiftCalculator(const JoukowskyTransform& jt, const VelocityCalculator& vc,
                               double rho_fluid, int n)
    : jt(jt), vc(vc), rho_fluid(rho_fluid), n(n) {}

std::complex<double> LiftCalculator::force(const std::vector<double>& psi) const {
    double r        = jt.circleRadius();
    double dtheta   = 2.0 * M_PI / n;
    double rho_eval = r * (1.0 + kDelta);  // just outside the boundary

    std::complex<double> F = 0.0;
    for (int i = 0; i < n; ++i) {
        double theta = i * dtheta;
        std::complex<double> e_ith = std::exp(std::complex<double>(0.0, theta));

        // velocity evalued just outside the boundary; close enough that it
	// approximates surface velocity but far enough that poisson riemann
	// sum can be taken
        std::complex<double> z_eval = jt.getB() + rho_eval * e_ith;
        auto q   = vc.velocity(z_eval, psi);
        double V_sq = std::norm(q);

        //   n_hat * ds = f'(z) * e^{i*theta} * r * dtheta
        // (the |f'(z)| in the unit normal and in ds cancel, so we use the
        //  geometric quantities at the actual boundary, not the offset point)
        std::complex<double> z_bnd  = jt.getB() + r * e_ith;
        std::complex<double> n_ds   = jt.dz2w(z_bnd) * e_ith * r * dtheta;

        F += 0.5 * rho_fluid * V_sq * n_ds; //bernoulli
    }
    return F;
}

double LiftCalculator::lift(const std::vector<double>& psi) const {
    return force(psi).imag();
}

double LiftCalculator::horizontalForce(const std::vector<double>& psi) const {
    return force(psi).real();
}
