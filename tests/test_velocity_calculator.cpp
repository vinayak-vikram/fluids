#include "../velocity_calculator.h"
#include "../boundary_builder.h"
#include "../poisson.h"
#include "../joukowsky.h"
#include <cmath>
#include <cstdio>
#include <complex>

static int failures = 0;

static void check(const char* name, bool ok) {
    if (ok) std::printf("  PASS  %s\n", name);
    else   { std::printf("  FAIL  %s\n", name); ++failures; }
}

static void check_near(const char* name, double got, double expected, double tol = 5e-3) {
    bool ok = std::abs(got - expected) < tol;
    if (!ok) std::printf("        got=%.6f  expected=%.6f\n", got, expected);
    check(name, ok);
}

// --- flat plate analytical case ---
//
// b=0, a=1, hinge=0 → generating circle = unit circle.
// w = z + 1/z, so on |z|=1: w = 2cosθ (flat plate on real axis).
//
// Exterior stream function: Ψ(ρ,θ) = ω(1 − cos(2θ)/ρ²)
//
// At z=2  (ρ=2, θ=0):  q_z = iω/4,  f'(2)=3/4  → q_w = iω/3
// At z=2i (ρ=2, θ=π/2):  q_z = -ω/4,  f'(2i)=5/4  → q_w = -ω/5

static const int    N     = 2000;
static const double OMEGA = 2.0;

static void test_flat_plate_velocity_z2() {
    JoukowskyTransform   jt(1.0, {0.0, 0.0});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, N);
    PoissonSolver        ps(jt.circleRadius());
    VelocityCalculator   vc(ps, jt);

    auto psi = bb.build(OMEGA);
    auto q   = vc.velocity({2.0, 0.0}, psi);
    check_near("flat plate z=2: Re(q) = 0",       q.real(),  0.0);
    check_near("flat plate z=2: Im(q) = omega/3", q.imag(),  OMEGA / 3.0);
}

static void test_flat_plate_velocity_z2i() {
    JoukowskyTransform   jt(1.0, {0.0, 0.0});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, N);
    PoissonSolver        ps(jt.circleRadius());
    VelocityCalculator   vc(ps, jt);

    auto psi = bb.build(OMEGA);
    auto q   = vc.velocity({0.0, 2.0}, psi);
    check_near("flat plate z=2i: Re(q) = -omega/5", q.real(), -OMEGA / 5.0);
    check_near("flat plate z=2i: Im(q) = 0",        q.imag(),  0.0);
}

static void test_velocity_scales_with_omega() {
    JoukowskyTransform   jt(1.0, {-0.1, 0.05});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, N);
    PoissonSolver        ps(jt.circleRadius());
    VelocityCalculator   vc(ps, jt);

    std::complex<double> z_test = {3.0, 1.0};
    auto q1 = vc.velocity(z_test, bb.build(1.0));
    auto q3 = vc.velocity(z_test, bb.build(3.0));

    check_near("velocity linear in omega (Re)", q3.real(), 3.0 * q1.real(), 1e-10);
    check_near("velocity linear in omega (Im)", q3.imag(), 3.0 * q1.imag(), 1e-10);
}

static void test_velocity_decays_far_from_wing() {
    JoukowskyTransform   jt(1.0, {-0.1, 0.0});
    BoundaryValueBuilder bb(jt, {0.0, 0.0}, N);
    PoissonSolver        ps(jt.circleRadius());
    VelocityCalculator   vc(ps, jt);

    auto psi    = bb.build(1.0);
    auto q_near = vc.velocity({2.0,  0.0}, psi);
    auto q_far  = vc.velocity({20.0, 0.0}, psi);

    check("|q| decays: far < near", std::abs(q_far) < std::abs(q_near));
    check("|q| far is small",       std::abs(q_far) < 1e-2);
}

int main() {
    std::printf("=== VelocityCalculator ===\n");
    test_flat_plate_velocity_z2();
    test_flat_plate_velocity_z2i();
    test_velocity_scales_with_omega();
    test_velocity_decays_far_from_wing();

    std::printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
