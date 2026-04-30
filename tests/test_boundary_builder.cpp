#include "../boundary_builder.h"
#include "../joukowsky.h"
#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;

static void check(const char* name, bool ok) {
    if (ok) std::printf("  PASS  %s\n", name);
    else   { std::printf("  FAIL  %s\n", name); ++failures; }
}

static void check_near(const char* name, double got, double expected, double tol = 1e-3) {
    bool ok = std::abs(got - expected) < tol;
    if (!ok) std::printf("        got=%.6f expected=%.6f\n", got, expected);
    check(name, ok);
}

// --- flat plate analytical case ---
//
// b=0, a=1, hinge=0 → generating circle is the unit circle.
// w = f(z) = z + 1/z, so on |z|=1: w = 2cosθ (real segment [-2,2]).
// f'(z) = 1 - 1/z²
// V_wing = iω(2cosθ - 0) = 2iω cosθ
//
// g(θ) = Re(V_wing · conj(f'(z)) · e^{-iθ})
//       = 2ω sin(2θ)                    ← zero mean: no correction needed
//
// ψ(θ) = r · ∫₀^θ 2ω sin(2φ) dφ = ω(1 - cos 2θ)
//
// Key values: ψ(π/2) = 2ω,  ψ(π) = 0,  ψ(3π/2) = 2ω

static void test_flat_plate_known_values() {
    JoukowskyTransform jt(1.0, {0.0, 0.0});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, 2000);
    double omega = 3.0;
    auto psi = bb.build(omega);
    int n = (int)psi.size();

    check_near("flat plate psi(0)    = 0",    psi[0],      0.0,     1e-10);
    check_near("flat plate psi(pi/2) = 2w",   psi[n / 4],  2*omega, 1e-2);
    check_near("flat plate psi(pi)   = 0",    psi[n / 2],  0.0,     1e-2);
    check_near("flat plate psi(3pi/2)= 2w",   psi[3*n/4],  2*omega, 1e-2);
}

static void test_zero_omega_gives_zero_psi() {
    JoukowskyTransform jt(1.0, {-0.1, 0.05});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, 500);
    auto psi = bb.build(0.0);
    bool all_zero = true;
    for (double v : psi) if (std::abs(v) > 1e-12) all_zero = false;
    check("zero omega → zero psi everywhere", all_zero);
}

static void test_output_size() {
    JoukowskyTransform jt(1.0, {-0.1, 0.0});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, 300);
    check("output size equals n", (int)bb.build(1.0).size() == 300);
}

static void test_linearity_in_omega() {
    JoukowskyTransform jt(1.5, {-0.1, 0.05});
    BoundaryValueBuilder bb(jt, {-0.5, 0.0}, 400);
    auto psi1 = bb.build(1.0);
    auto psi2 = bb.build(3.0);
    bool linear = true;
    for (int i = 0; i < 400; ++i)
        if (std::abs(psi2[i] - 3.0 * psi1[i]) > 1e-10) linear = false;
    check("psi scales linearly with omega", linear);
}

static void test_periodicity() {
    JoukowskyTransform jt(1.0, {-0.08, 0.04});
    int n = 800;
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, n);
    double r = bb.circleRadius();
    auto psi = bb.build(2.0);
    check_near("periodicity: psi[0] = 0",      psi[0],   0.0, 1e-10);
    check("periodicity: psi bounded by r*2pi", std::abs(psi[n/4]) < r * 2.0 * M_PI);
}

static void test_circle_radius() {
    JoukowskyTransform jt(2.0, {-0.3, 0.1});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, 100);
    double expected = jt.circleRadius();
    check_near("circleRadius() delegates to jt", bb.circleRadius(), expected, 1e-12);
}

int main() {
    std::printf("=== BoundaryValueBuilder ===\n");
    test_zero_omega_gives_zero_psi();
    test_output_size();
    test_linearity_in_omega();
    test_circle_radius();
    test_periodicity();
    test_flat_plate_known_values();

    std::printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
