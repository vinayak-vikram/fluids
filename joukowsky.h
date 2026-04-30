#ifndef JOUKOWSKY_H
#define JOUKOWSKY_H

#include <complex>

class JoukowskyTransform {
public:
    JoukowskyTransform(double a, std::complex<double> b);
    std::complex<double> z2w(std::complex<double> z) const;
    std::complex<double> w2z(std::complex<double> w) const; // returns two solns
    std::complex<double> dz2w(std::complex<double> z) const; // derivative dw/dz

    double               getA() const { return a; }
    std::complex<double> getB() const { return b; }
    double               circleRadius() const { return std::abs(std::complex<double>(a, 0.0) - b); }

private:
    double a;
    std::complex<double> b;
};

#endif
