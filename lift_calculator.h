#ifndef LIFT_CALCULATOR_H
#define LIFT_CALCULATOR_H

#include "joukowsky.h"
#include "velocity_calculator.h"
#include <complex>
#include <vector>

class LiftCalculator {
public:
    LiftCalculator(const JoukowskyTransform& jt, const VelocityCalculator& vc,
                   double rho_fluid, int n);
    double lift(const std::vector<double>& psi) const;
    double horizontalForce(const std::vector<double>& psi) const;

private:
    std::complex<double> force(const std::vector<double>& psi) const;

    const JoukowskyTransform& jt;
    const VelocityCalculator& vc;
    double rho_fluid;
    int    n;

    static constexpr double kDelta = 0.02;
};

#endif
