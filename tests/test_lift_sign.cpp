#include "../lift_calculator.h"
#include "../velocity_calculator.h"
#include "../boundary_builder.h"
#include "../poisson.h"
#include "../joukowsky.h"
#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>
#include <algorithm>

// Full pipeline in one struct.
struct Pipeline {
    JoukowskyTransform   jt;
    PoissonSolver        ps;
    VelocityCalculator   vc;
    BoundaryValueBuilder bb;
    LiftCalculator       lc;

    Pipeline(double a, std::complex<double> b, std::complex<double> hinge,
             int n, double rho)
        : jt(a, b), ps(jt.circleRadius()), vc(ps, jt)
        , bb(jt, hinge, n), lc(jt, vc, rho, n) {}
};

// ── helpers ──────────────────────────────────────────────────────────────────
static int failures = 0;
static void check(const char* name, bool ok) {
    std::printf("  %s  %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) ++failures;
}

// ── find leading edge (leftmost w-plane point on the generating circle) ──────
static std::complex<double> leading_edge(const JoukowskyTransform& jt) {
    double r0 = jt.circleRadius();
    std::complex<double> best = jt.z2w(jt.getB() + r0);
    for (int i = 0; i < 1000; i++) {
        double theta = 2*M_PI*i/1000.0;
        auto z = jt.getB() + r0 * std::exp(std::complex<double>(0, theta));
        auto w = jt.z2w(z);
        if (w.real() < best.real()) best = w;
    }
    return best;
}

static std::complex<double> trailing_edge(const JoukowskyTransform& jt) {
    // Critical point z = a
    return jt.z2w(std::complex<double>(jt.getA(), 0.0));
}

// ── test 1: print raw lift values for different hinges ────────────────────────
static void test_lift_sign_vs_hinge() {
    const double a = 1.0;
    const std::complex<double> b{-0.1, 0.08};
    const int N = 500;
    const double omega = 2.0;
    const double rho = 1.225;

    JoukowskyTransform jt(a, b);
    auto h_lead = leading_edge(jt);
    auto h_trail = trailing_edge(jt);
    auto h_ctr = std::complex<double>(0.0, 0.0);

    std::printf("\n=== Lift sign vs hinge (cambered airfoil B={%.2f,%.2f}) ===\n",
        b.real(), b.imag());
    std::printf("  Leading edge hinge:  w = (%.4f, %.4f)\n", h_lead.real(), h_lead.imag());
    std::printf("  Trailing edge hinge: w = (%.4f, %.4f)\n", h_trail.real(), h_trail.imag());
    std::printf("  Omega = %.1f rad/s\n\n", omega);

    auto print_lift = [&](const char* label, std::complex<double> hinge) {
        Pipeline p(a, b, hinge, N, rho);
        auto psi = p.bb.build(omega);
        double L = p.lc.lift(psi);
        double H = p.lc.horizontalForce(psi);
        std::printf("  %-24s L = %+10.4f  H = %+10.4f  |F| = %.4f\n",
            label, L, H, std::hypot(L, H));
    };

    print_lift("leading edge hinge:", h_lead);
    print_lift("trailing edge hinge:", h_trail);
    print_lift("center hinge (0,0):", h_ctr);
    print_lift("quarter-chord hinge:", (3.0*h_lead + h_trail) / 4.0);
}

// ── test 2: verify normal vectors point outward (Im(n_ds) > 0 on top arc) ────
static void test_normals_outward() {
    const double a = 1.0;
    const std::complex<double> b{0.0, 0.0};  // flat plate for clarity
    JoukowskyTransform jt(a, b);

    std::printf("\n=== Normal direction check (flat plate B=0) ===\n");

    double r = jt.circleRadius();
    // At theta = pi/2 (top of circle, upper surface of flat plate)
    double theta = M_PI / 2.0;
    auto z = jt.getB() + r * std::exp(std::complex<double>(0, theta));
    auto fprime = jt.dz2w(z);
    auto e_ith = std::exp(std::complex<double>(0, theta));
    auto n_ds = fprime * e_ith * r;  // without dtheta

    std::printf("  At theta=pi/2 (top):    n_ds = (%+.4f, %+.4f)  Im>0? %s\n",
        n_ds.real(), n_ds.imag(), n_ds.imag() > 0 ? "YES (outward=up)" : "NO");

    // At theta = 3*pi/2 (bottom)
    theta = 3.0 * M_PI / 2.0;
    z = jt.getB() + r * std::exp(std::complex<double>(0, theta));
    fprime = jt.dz2w(z);
    e_ith = std::exp(std::complex<double>(0, theta));
    n_ds = fprime * e_ith * r;

    std::printf("  At theta=3pi/2 (bot):   n_ds = (%+.4f, %+.4f)  Im<0? %s\n",
        n_ds.real(), n_ds.imag(), n_ds.imag() < 0 ? "YES (outward=down)" : "NO");

    check("outward normal is up on top surface",
          [&]{ double th=M_PI/2; auto z2=jt.getB()+r*std::exp(std::complex<double>(0,th));
               return (jt.dz2w(z2)*std::exp(std::complex<double>(0,th))*r).imag() > 0; }());
    check("outward normal is down on bottom surface",
          [&]{ double th=3*M_PI/2; auto z2=jt.getB()+r*std::exp(std::complex<double>(0,th));
               return (jt.dz2w(z2)*std::exp(std::complex<double>(0,th))*r).imag() < 0; }());
}

// ── test 3: velocity on upper vs lower surface (leading edge hinge) ───────────
// For positive lift: fluid must be faster on the upper surface.
static void test_velocity_asymmetry() {
    const double a = 1.0;
    const std::complex<double> b{-0.1, 0.08};
    const int N = 800;
    const double omega = 2.0;

    JoukowskyTransform   jt(a, b);
    PoissonSolver        ps(jt.circleRadius());
    VelocityCalculator   vc(ps, jt);

    auto h_lead  = leading_edge(jt);
    auto h_trail = trailing_edge(jt);

    auto sample_surfaces = [&](std::complex<double> hinge, const char* label) {
        BoundaryValueBuilder bb(jt, hinge, N);
        auto psi = bb.build(omega);

        // Sample 8 evenly spaced points on upper surface (theta in (0,pi))
        // and 8 on lower surface (theta in (pi,2pi)), slightly outside boundary.
        double r0 = jt.circleRadius();
        double rho_eval = r0 * 1.03;

        double upper_sum = 0, lower_sum = 0;
        int ns = 8;
        for (int i = 0; i < ns; i++) {
            // upper: theta in (0, pi)
            double th_up = (i + 0.5) * M_PI / ns;
            auto z_up = jt.getB() + rho_eval * std::exp(std::complex<double>(0, th_up));
            double v_up = std::abs(vc.velocity(z_up, psi));

            // lower: theta in (pi, 2pi) i.e. (pi, pi + (i+0.5)*pi/ns)
            double th_lo = M_PI + (i + 0.5) * M_PI / ns;
            auto z_lo = jt.getB() + rho_eval * std::exp(std::complex<double>(0, th_lo));
            double v_lo = std::abs(vc.velocity(z_lo, psi));

            upper_sum += v_up;
            lower_sum += v_lo;
        }
        double avg_upper = upper_sum / ns;
        double avg_lower = lower_sum / ns;
        std::printf("  %-28s avg|V|_upper=%6.4f  avg|V|_lower=%6.4f  upper>lower? %s\n",
            label, avg_upper, avg_lower,
            avg_upper > avg_lower ? "YES → positive lift expected" : "NO → negative lift expected");
    };

    std::printf("\n=== Velocity asymmetry (upper vs lower surface) ===\n");
    sample_surfaces(h_lead,  "leading edge hinge:");
    sample_surfaces(h_trail, "trailing edge hinge:");
}

// ── test 4: sign assertion — cambered airfoil with leading edge hinge ─────────
// Physical expectation: with ω>0 (trailing edge swings UP), the induced flow
// is faster on the upper (cambered) surface → lower Bernoulli pressure above
// → net upward force → positive lift.
static void test_leading_edge_hinge_positive_lift() {
    const double a = 1.0;
    const std::complex<double> b{-0.1, 0.08};
    const int N = 1000;
    const double omega = 2.0;
    const double rho = 1.225;

    JoukowskyTransform jt(a, b);
    auto h_lead = leading_edge(jt);
    Pipeline p(a, b, h_lead, N, rho);
    auto psi = p.bb.build(omega);
    double L = p.lc.lift(psi);

    std::printf("\n  Leading edge hinge lift at omega=%.1f: L = %+.6f\n", omega, L);
    check("leading edge hinge → positive lift", L > 0.0);
}

// ── test 5: Bernoulli pressure sign check ────────────────────────────────────
// Explicitly verify: area of higher |V|^2 on upper surface gives positive Im(F).
static void test_bernoulli_sign_manual() {
    // Construct a fake psi that has |V| larger on upper surface.
    // Use a simple harmonic: psi(theta) = sin(theta) → more flow on upper half.
    // Then manually compute the Bernoulli integral and check sign.
    //
    // For psi = A*sin(theta) on a circle of radius r:
    // The exterior solution: Psi(rho,theta) = A*(r/rho)*sin(theta)
    // dPsi/drho = -A*r/rho^2 * sin(theta)
    // dPsi/dtheta = A*r/rho * cos(theta)
    // q_z = e^{-ith}*(i*dPsi/drho + dPsi/dtheta/rho)
    //      = e^{-ith}*(i*(-A*r/rho^2)*sin(th) + A*r/rho^2*cos(th))
    //      = (A*r/rho^2)*e^{-ith}*(cos(th) - i*sin(th))
    //      = (A*r/rho^2)*e^{-ith}*e^{-ith} = (A*r/rho^2)*e^{-2ith}
    //
    // For flat plate (B=0, a=1): f'(z) = 1 - 1/z^2
    // At rho=r=1 (boundary), z = e^{ith}: q_z = A*e^{-2ith}
    // q_w = q_z / f'(z) = A*e^{-2ith} / (1 - e^{-2ith})
    // This diverges at theta=0 (trailing edge) — expected for potential flow.
    //
    // Instead of computing analytically, verify numerically using PoissonSolver.
    const double a = 1.0;
    const std::complex<double> b{0.0, 0.0};  // flat plate
    const int N = 2000;
    const double r0 = std::abs(std::complex<double>(a,0.0) - b);  // = 1

    // Build psi = sin(theta): upper half has positive psi, lower has negative
    std::vector<double> psi(N);
    for (int i = 0; i < N; i++) {
        double theta = 2*M_PI*i/N;
        psi[i] = std::sin(theta);  // positive on upper half
    }

    JoukowskyTransform   jt(a, b);
    PoissonSolver        ps(r0);
    VelocityCalculator   vc(ps, jt);
    LiftCalculator       lc(jt, vc, 1.0, N);  // rho=1 for simplicity

    double L = lc.lift(psi);
    double H = lc.horizontalForce(psi);
    std::printf("\n=== Bernoulli sign manual check (flat plate, psi=sin(theta)) ===\n");
    std::printf("  L = %+.6f  H = %+.6f\n", L, H);
    std::printf("  psi=sin(theta) concentrates flow on upper surface (theta in 0..pi)\n");
    // The normal on upper surface has Im(n_ds)>0 (upward), so if |V| is larger there,
    // Im(F) = L should be positive.
    check("sin(theta) boundary → positive lift (upward force)", L > 0.0);
}

int main() {
    std::printf("=== LiftSign diagnostic ===\n");
    test_normals_outward();
    test_bernoulli_sign_manual();
    test_velocity_asymmetry();
    test_lift_sign_vs_hinge();
    test_leading_edge_hinge_positive_lift();

    std::printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
