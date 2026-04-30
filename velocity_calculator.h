#ifndef VELOCITY_CALCULATOR_H
#define VELOCITY_CALCULATOR_H

#include "poisson.h"
#include "joukowsky.h"
#include <complex>
#include <vector>

class VelocityCalculator {
public:
    VelocityCalculator(const PoissonSolver& ps, const JoukowskyTransform& jt);
    std::complex<double> velocity(std::complex<double> z,
                                  const std::vector<double>& psi) const;

private:
    const PoissonSolver&      ps;
    const JoukowskyTransform& jt;

    static constexpr double kH = 1e-4;
};

#endif
