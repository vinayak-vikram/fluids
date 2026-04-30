#include "poisson.h"

PoissonSolver::PoissonSolver(double radius) : r0(radius) {}

double PoissonSolver::calculateKernel(double r, double theta, double phi) const {
    double r0_sq = r0 * r0;
    double r_sq = r * r;
    double diff = phi - theta;
    
    double num = r0_sq - r_sq;
    double den = r0_sq - (2 * r0 * r * std::cos(diff)) + r_sq;
   
    // point on boundary causes div. by 0, fallback by kernel properties
    if (std::abs(den) < 1e-9) return 0.0;

    return num / den;
}

double PoissonSolver::calculateExteriorKernel(double r, double theta, double phi) const {
    double r0_sq = r0 * r0;
    double r_sq = r * r;
    double diff = phi - theta;

    double num = r_sq - r0_sq;
    double den = r0_sq - (2 * r0 * r * std::cos(diff)) + r_sq;

    if (std::abs(den) < 1e-9) return 0.0;

    return num / den;
}

double PoissonSolver::solveExternal(double r, double theta, const std::vector<double>& boundaryValues) const {
    int n = boundaryValues.size();
    double dphi = (2.0 * M_PI) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double phi = i * dphi;
        sum += calculateExteriorKernel(r, theta, phi) * boundaryValues[i] * dphi;
    }
    return sum / (2.0 * M_PI);
}

double PoissonSolver::solveInternal(double r, double theta, const std::vector<double>& boundaryValues) const {
    int n = boundaryValues.size();
    double dphi = (2.0 * M_PI) / n;
    double sum = 0.0;
    //riemann sum impl, should be computationally pretty doable for smallish n?
    for (int i = 0; i < n; ++i) {
        double phi = i * dphi; //quant
        double dP = calculateKernel(r, theta, phi);
	sum += dP * boundaryValues[i] * dphi;
    }
    return sum / (2.0 * M_PI);
}
