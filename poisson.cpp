#include "poisson.h"
#include <numeric>
#include <numbers> // pi

PoissonSolver::PoissonSolver(double radius) : r0(radius) {}

double PoissonSolver::calculateKernel(double r, double theta, double phi) const {
    double r0_sq = r0 * r0;
    double r_sq = r * r;
    double diff = phi - theta;
    
    double num = r0_sq - r_sq;
    double den = r0_sq - (2 * r0 * r * std::cos(diff)) + r_sq;
   
    // point on boundary causes div. by 0, fallback by kernel properties
    if (std::abs(denominator) < 1e-9) return 0.0; 

    return numerator / denominator;
}

double PoissonSolver::solveInternal(double r, double theta, const std::vector<double>& boundaryValues) const {
    int n = boundaryValues.size();
    double dphi = (2.0 * std::numbers::pi) / n;
    double sum = 0.0;
    //riemann sum impl, should be computationally pretty doable for smallish n?
    for (int i = 0; i < n; ++i) {
        double phi = i * dphi; //curr angle on the boundary
        double dP = calculateKernel(r, theta, phi);
	sum += dP * boundaryValues[i] * dphi;
    }
    return sum / (2.0 * std::numbers::pi);
}
