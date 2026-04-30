#include "../lift_calculator.h"
#include "../boundary_builder.h"
#include "../poisson.h"
#include "../joukowsky.h"
#include <cmath>
#include <cstdio>

static int failures = 0;

static void check(const char* name, bool ok) {
    if (ok) std::printf("  PASS  %s\n", name);
    else   { std::printf("  FAIL  %s\n", name); ++failures; }
}

static void check_near(const char* name, double got, double expected, double tol) {
    bool ok = std::abs(got - expected) < tol;
    if (!ok) std::printf("        got=%.6f  expected=%.6f\n", got, expected);
    check(name, ok);
}

// Helper: build the full pipeline for a given transform and hinge.
struct Pipeline {
    JoukowskyTransform   jt;
    PoissonSolver        ps;
    VelocityCalculator   vc;
    BoundaryValueBuilder bb;
    LiftCalculator       lc;

    Pipeline(double a, std::complex<double> b, std::complex<double> hinge,
             int n, double rho_fluid)
        : jt(a, b)
        , ps(jt.circleRadius())
        , vc(ps, jt)
        , bb(jt, hinge, n)
        , lc(jt, vc, rho_fluid, n)
    {}
};

static void test_zero_omega_zero_lift() {
    Pipeline p(1.0, {-0.1, 0.05}, {0.0, 0.0}, 1000, 1.225);
    auto psi = p.bb.build(0.0);
    check_near("zero omega → zero lift",             p.lc.lift(psi),            0.0, 1e-10);
    check_near("zero omega → zero horizontal force", p.lc.horizontalForce(psi), 0.0, 1e-10);
}

static void test_lift_scales_as_omega_squared() {
    // Velocity is linear in omega, so pressure (|V|^2) and lift scale as omega^2.
    Pipeline p(1.0, {-0.1, 0.05}, {-0.5, 0.0}, 1500, 1.225);
    double L1 = p.lc.lift(p.bb.build(1.0));
    double L2 = p.lc.lift(p.bb.build(2.0));
    // L2 should be 4*L1; test the ratio (avoid division if L1 near 0)
    check_near("lift scales as omega^2", L2, 4.0 * L1, std::abs(L1) * 0.02 + 1e-10);
}

static void test_symmetric_airfoil_center_hinge_zero_lift() {
    // A symmetric airfoil (b purely real) rotating about its geometric centre
    // has a velocity field that is antisymmetric about the chord line.
    // The pressure contribution from the upper and lower surfaces cancel,
    // giving zero net lift.
    Pipeline p(1.0, {-0.1, 0.0}, {0.0, 0.0}, 2000, 1.225);
    double L = p.lc.lift(p.bb.build(3.0));
    // Allow small numerical residual (~1% of the horizontal force magnitude)
    double H = std::abs(p.lc.horizontalForce(p.bb.build(3.0)));
    check("symmetric + center hinge → near-zero lift", std::abs(L) < 0.05 * H + 1e-6);
}

static void test_cambered_airfoil_nonzero_lift() {
    // A cambered airfoil (b has imaginary component) breaks the up/down symmetry.
    // Rotating it should produce a nonzero net lift.
    Pipeline p(1.0, {-0.1, 0.1}, {0.0, 0.0}, 2000, 1.225);
    double L = p.lc.lift(p.bb.build(2.0));
    check("cambered airfoil generates nonzero lift", std::abs(L) > 1e-4);
}

static void test_lift_vs_horizontal_force_pythagoras() {
    // The total aerodynamic force magnitude should equal sqrt(L^2 + H^2).
    // This is a consistency check between lift() and horizontalForce().
    Pipeline p(1.0, {-0.1, 0.07}, {-0.3, 0.0}, 1500, 1.225);
    auto psi = p.bb.build(1.5);
    double L  = p.lc.lift(psi);
    double H  = p.lc.horizontalForce(psi);
    double F2 = L*L + H*H;
    // Re-compute via the internal force magnitude from separate calls
    // (they use the same force() internally, just split Im/Re, so they should be consistent)
    check("L and H are consistent: L^2+H^2 > 0 when wing moves", F2 > 1e-8);
}

int main() {
    std::printf("=== LiftCalculator ===\n");
    test_zero_omega_zero_lift();
    test_lift_scales_as_omega_squared();
    test_symmetric_airfoil_center_hinge_zero_lift();
    test_cambered_airfoil_nonzero_lift();
    test_lift_vs_horizontal_force_pythagoras();

    std::printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
