#ifndef BOUNDARY_BUILDER_H
#define BOUNDARY_BUILDER_H

#include "joukowsky.h"
#include <complex>
#include <vector>

class BoundaryValueBuilder {
public:
    BoundaryValueBuilder(const JoukowskyTransform& jt,
                         std::complex<double> hinge, int n);
    // stream-function values
    std::vector<double> build(double omega) const;

    double circleRadius() const { return jt.circleRadius(); }

private:
    const JoukowskyTransform& jt;
    std::complex<double>      hinge;
    int                       n;
};

#endif
