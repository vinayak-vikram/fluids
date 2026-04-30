#include "joukowsky.h"

JoukowskyTransform::JoukowskyTransform(double a, std::complex<double> b) : a(a), b(b) {}

std::complex<double> JoukowskyTransform::z2w(std::complex<double> z) const {
    return z + (a * a) / z;
}

// dw/dz = 1 - a^2/z^2
std::complex<double> JoukowskyTransform::dz2w(std::complex<double> z) const {
    return 1.0 - (a * a) / (z * z);
}

// solve z^2-wz+a^2=0
// returns the exterior root, subject to change but should suffice for now?
std::complex<double> JoukowskyTransform::w2z(std::complex<double> w) const {
    std::complex<double> disc = std::sqrt(w * w - 4.0 * a * a);
    std::complex<double> z1 = (w + disc) / 2.0;
    std::complex<double> z2 = (w - disc) / 2.0;
    return std::abs(z1) >= std::abs(z2) ? z1 : z2;
}
