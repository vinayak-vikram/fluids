#include "boundary_builder.h"
#include <cmath>
#include <numeric>

BoundaryValueBuilder::BoundaryValueBuilder(const JoukowskyTransform& jt,
                                           double a, std::complex<double> b,
                                           std::complex<double> hinge, int n)
    : jt(jt), a(a), b(b), hinge(hinge), n(n) {}

double BoundaryValueBuilder::circleRadius() const {
    return std::abs(std::complex<double>(a, 0.0) - b);
}

std::vector<double> BoundaryValueBuilder::build(double omega) const {
    double r      = circleRadius();
    double dtheta = 2.0 * M_PI / n;

    // compute neumann data g(theta_1) in the z-plane, see textbook
    std::vector<double> g(n);
    for (int i = 0; i < n; ++i) {
        double theta = i * dtheta;
        std::complex<double> z      = b + r * std::exp(std::complex<double>(0.0, theta));
        std::complex<double> w      = jt.z2w(z);
        std::complex<double> fprime = 1.0 - (a * a) / (z * z);
        std::complex<double> V_wing = std::complex<double>(0.0, omega) * (w - hinge);
        std::complex<double> e_neg  = std::exp(std::complex<double>(0.0, -theta));
        g[i] = (V_wing * std::conj(fprime) * e_neg).real();
    }

    // enforce zero net flux
    // analytically exact; subtract mean for robustness
    double mean_g = std::accumulate(g.begin(), g.end(), 0.0) / n;
    for (double& gi : g) gi -= mean_g;

    // integrate to get stream-function boundary values
    // just cre lmao
    std::vector<double> psi(n);
    psi[0] = 0.0;
    for (int i = 1; i < n; ++i)
        psi[i] = psi[i - 1] + r * g[i - 1] * dtheta;

    return psi;
}
