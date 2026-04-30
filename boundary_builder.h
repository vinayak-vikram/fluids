#ifndef BOUNDARY_BUILDER_H
#define BOUNDARY_BUILDER_H

#include "joukowsky.h"
#include <complex>
#include <vector>

class BoundaryValueBuilder {
public:
    BoundaryValueBuilder(const JoukowskyTransform& jt,
                         double a, std::complex<double> b,
                         std::complex<double> hinge, int n);
    std::vector<double> build(double omega) const;
    double circleRadius() const;

private:
    const JoukowskyTransform& jt;
    double a;
    std::complex<double> b;
    std::complex<double> hinge;
    int n;
};

#endif
