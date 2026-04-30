#include "velocity_calculator.h"
#include <cmath>

VelocityCalculator::VelocityCalculator(const PoissonSolver& ps,
                                       const JoukowskyTransform& jt)
    : ps(ps), jt(jt) {}

std::complex<double> VelocityCalculator::velocity(std::complex<double> z,
                                                   const std::vector<double>& psi) const {
    //rel to b
    std::complex<double> zeta = z - jt.getB();
    double rho   = std::abs(zeta);
    double theta = std::arg(zeta);

    //solve via central differences
    double dPsi_drho = (ps.solveExternal(rho + kH, theta, psi)
                      - ps.solveExternal(rho - kH, theta, psi)) / (2.0 * kH);

    double dPsi_dtheta = (ps.solveExternal(rho, theta + kH, psi)
                        - ps.solveExternal(rho, theta - kH, psi)) / (2.0 * kH);

    //use polar form of cre to get complex velocity (z-plane)
    std::complex<double> e_neg_itheta = std::exp(std::complex<double>(0.0, -theta));
    std::complex<double> q_z = e_neg_itheta
                              * (std::complex<double>(0.0, 1.0) * dPsi_drho
                                 + dPsi_dtheta / rho);

    // pull through conformal factor
    return q_z / jt.dz2w(z);
}
