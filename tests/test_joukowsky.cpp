#include "../joukowsky.h"
#include <cmath>
#include <cstdio>
#include <complex>

static int failures = 0;

static void check(const char* name, bool ok) {
    if (ok) std::printf("  PASS  %s\n", name);
    else   { std::printf("  FAIL  %s\n", name); ++failures; }
}

static void check_near(const char* name, std::complex<double> got, std::complex<double> expected, double tol = 1e-9) {
    check(name, std::abs(got - expected) < tol);
}

// z2w(z) = z + a²/z
// For a=1, z=2: w = 2 + 1/2 = 2.5
// For a=1, z=i: w = i + 1/i = i - i = 0  (flat plate collapses imaginary axis)
static void test_z2w_known() {
    JoukowskyTransform jt(1.0, {0.0, 0.0});
    check_near("z2w a=1 z=2",      jt.z2w({2.0, 0.0}), {2.5, 0.0});
    check_near("z2w a=1 z=i",      jt.z2w({0.0, 1.0}), {0.0, 0.0});
    check_near("z2w a=1 z=-i",     jt.z2w({0.0,-1.0}), {0.0, 0.0});
    // On the critical circle |z|=a the map folds; z=a maps to w=2a
    check_near("z2w critical pt",  jt.z2w({1.0, 0.0}), {2.0, 0.0});
}

// w2z should be the inverse of z2w for |z| > a
static void test_roundtrip_exterior() {
    JoukowskyTransform jt(1.5, {-0.1, 0.05});
    // Sample a few points clearly outside the critical circle
    std::complex<double> pts[] = {{3.0, 1.0}, {2.0, -0.5}, {4.0, 2.0}, {1.8, 0.3}};
    for (auto& z : pts) {
        std::complex<double> w = jt.z2w(z);
        std::complex<double> zr = jt.w2z(w);
        check_near("roundtrip z→w→z", zr, z);
    }
}

// The two roots of the inverse multiply to a² (Vieta's)
static void test_w2z_vieta() {
    JoukowskyTransform jt(2.0, {0.0, 0.0});
    std::complex<double> w = {5.0, 1.0};
    std::complex<double> z_ext = jt.w2z(w);
    std::complex<double> z_int = (2.0 * 2.0) / z_ext;  // other root = a²/z_ext
    // Both roots should satisfy the map
    check_near("w2z exterior root satisfies map", jt.z2w(z_ext), w);
    check_near("w2z interior root satisfies map", jt.z2w(z_int), w);
}

// w2z returns exterior root: |z| >= a
static void test_w2z_returns_exterior() {
    JoukowskyTransform jt(1.0, {0.0, 0.0});
    std::complex<double> ws[] = {{3.0, 0.5}, {-2.0, 1.0}, {0.0, 4.0}};
    for (auto& w : ws)
        check("w2z |z| >= a", std::abs(jt.w2z(w)) >= 1.0 - 1e-12);
}

// Airfoil profile: circle |z - b| = |a - b| traced through z2w should be symmetric
// about the real axis when b is purely real (symmetric airfoil, no camber)
static void test_symmetric_airfoil_real_axis() {
    double a = 1.0;
    std::complex<double> b = {-0.1, 0.0};  // real offset → symmetric airfoil
    JoukowskyTransform jt(a, b);
    double r = std::abs(std::complex<double>(a, 0.0) - b);  // circle radius

    int n = 64;
    bool symmetric = true;
    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * i / n;
        std::complex<double> z_top = b + r * std::exp(std::complex<double>(0.0,  theta));
        std::complex<double> z_bot = b + r * std::exp(std::complex<double>(0.0, -theta));
        std::complex<double> w_top = jt.z2w(z_top);
        std::complex<double> w_bot = jt.z2w(z_bot);
        // Real parts match, imaginary parts are negated
        if (std::abs(w_top.real() - w_bot.real()) > 1e-10) symmetric = false;
        if (std::abs(w_top.imag() + w_bot.imag()) > 1e-10) symmetric = false;
    }
    check("symmetric airfoil (real b) is symmetric about real axis", symmetric);
}

// Cambered airfoil: imaginary b offset shifts the circle, making the airfoil
// asymmetric about the real axis — max_y + min_y != 0 (nonzero mean camber).
static void test_cambered_airfoil_is_asymmetric() {
    double a = 1.0;
    std::complex<double> b = {0.0, 0.05};
    JoukowskyTransform jt(a, b);
    double r = std::abs(std::complex<double>(a, 0.0) - b);

    double max_y = -1e9, min_y = 1e9;
    int n = 200;
    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * i / n;
        std::complex<double> z = b + r * std::exp(std::complex<double>(0.0, theta));
        double y = jt.z2w(z).imag();
        max_y = std::max(max_y, y);
        min_y = std::min(min_y, y);
    }
    // Mean camber height is nonzero (airfoil not symmetric about chord line)
    check("cambered airfoil has nonzero mean camber (max_y + min_y != 0)",
          std::abs(max_y + min_y) > 1e-3);
}

int main() {
    std::printf("=== Joukowsky transform ===\n");
    test_z2w_known();
    test_roundtrip_exterior();
    test_w2z_vieta();
    test_w2z_returns_exterior();
    test_symmetric_airfoil_real_axis();
    test_cambered_airfoil_is_asymmetric();

    std::printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
