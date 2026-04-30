#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <cmath>

class PoissonSolver {
public:
    PoissonSolver(double r0);
    double calculateKernel(double r, double theta, double phi) const;
    double calculateExteriorKernel(double r, double theta, double phi) const;
    // boundaryValues is a vector representing u(r0, phi) sampled around the circle
    double solveInternal(double r, double theta, const std::vector<double>& boundaryValues) const;
    double solveExternal(double r, double theta, const std::vector<double>& boundaryValues) const;
private:
    double r0;
};

#endif
