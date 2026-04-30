#ifndef JOUKOWSKY_H
#define JOUKOWSKY_H

#include <complex>

class JoukowskyTransform {
public:
    JoukowskyTransform(double a, std::complex<double> b);
    std::complex<double> z2w(std::complex<double> z) const;
    std::complex<double> w2z(std::complex<double> w) const; //returns two solns

private:
    double a;
    std::complex<double> b;
};

#endif
