#include "../poisson.h"
#include <cmath>
#include <cstdio>
#include <vector>

static int failures = 0;

static void check(const char* name, bool ok) {
    if (ok) {
        std::printf("  PASS  %s\n", name);
    } else {
        std::printf("  FAIL  %s\n", name);
        ++failures;
    }
}

static void check_near(const char* name, double got, double expected, double tol = 1e-6) {
    check(name, std::abs(got - expected) < tol);
}

// --- kernel tests ---

static void test_kernel_at_origin() {
    // K(0, θ, φ) = r0² / r0² = 1 for any angles
    PoissonSolver ps(2.0);
    check_near("kernel at origin (theta=0, phi=0)", ps.calculateKernel(0.0, 0.0, 0.0), 1.0);
    check_near("kernel at origin (theta=1, phi=3)", ps.calculateKernel(0.0, 1.0, 3.0), 1.0);
}

static void test_kernel_known_value() {
    // K(r, θ, φ) = (r0²-r²) / (r0²-2r0·r·cos(φ-θ)+r²)
    // With r0=2, r=1, θ=0, φ=0: (4-1)/(4-4+1) = 3/1 = 3
    PoissonSolver ps(2.0);
    check_near("kernel known value r=1 theta=phi=0", ps.calculateKernel(1.0, 0.0, 0.0), 3.0);

    // With r0=1, r=0.5, θ=0, φ=π: cos(π)=-1
    // num = 1-0.25 = 0.75, den = 1 - 2*0.5*(-1) + 0.25 = 2.25
    PoissonSolver ps2(1.0);
    check_near("kernel known value r=0.5 phi-theta=pi", ps2.calculateKernel(0.5, 0.0, M_PI), 0.75 / 2.25);
}

static void test_kernel_symmetry() {
    // K(r, θ, φ) depends only on φ-θ
    PoissonSolver ps(3.0);
    double k1 = ps.calculateKernel(1.0, 0.5, 1.2);
    double k2 = ps.calculateKernel(1.0, 0.5 + 1.0, 1.2 + 1.0);
    check_near("kernel angular translation symmetry", k1, k2);
}

static void test_kernel_positive_inside() {
    // For r < r0, kernel must be positive
    PoissonSolver ps(5.0);
    bool all_pos = true;
    for (int i = 0; i < 8; ++i) {
        double phi = i * M_PI / 4.0;
        if (ps.calculateKernel(3.0, 0.0, phi) <= 0.0) all_pos = false;
    }
    check("kernel positive for r < r0", all_pos);
}

static void test_kernel_boundary_singularity() {
    // r == r0 and theta == phi → denominator 0 → returns 0 (our fallback)
    PoissonSolver ps(2.0);
    check_near("kernel boundary singularity returns 0", ps.calculateKernel(2.0, 1.0, 1.0), 0.0);
}

// --- integral tests ---

static void test_integral_constant_boundary() {
    // Mean-value property: constant boundary c → solution is c everywhere inside
    PoissonSolver ps(1.0);
    int n = 1000;
    std::vector<double> boundary(n, 7.0);
    check_near("integral constant boundary r=0.5", ps.solveInternal(0.5, 0.3, boundary), 7.0, 1e-4);
    check_near("integral constant boundary r=0",   ps.solveInternal(0.0, 0.0, boundary), 7.0, 1e-4);
}

static void test_integral_at_origin_equals_mean() {
    // At r=0 kernel is 1 everywhere, so result = mean(boundary)
    PoissonSolver ps(1.0);
    int n = 500;
    std::vector<double> boundary(n);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        boundary[i] = std::sin(i * 2.0 * M_PI / n) + 2.0;
        sum += boundary[i];
    }
    double mean = sum / n;
    check_near("integral at origin equals mean of boundary", ps.solveInternal(0.0, 0.0, boundary), mean, 1e-10);
}

static void test_integral_theta_independence_at_origin() {
    // At r=0 the result must not depend on theta
    PoissonSolver ps(2.0);
    int n = 200;
    std::vector<double> boundary(n);
    for (int i = 0; i < n; ++i)
        boundary[i] = std::cos(i * 2.0 * M_PI / n) * 3.0 + 5.0;
    double v0 = ps.solveInternal(0.0, 0.0, boundary);
    double v1 = ps.solveInternal(0.0, 1.2, boundary);
    double v2 = ps.solveInternal(0.0, 2.7, boundary);
    check_near("integral at origin independent of theta (0 vs 1.2)", v0, v1, 1e-10);
    check_near("integral at origin independent of theta (0 vs 2.7)", v0, v2, 1e-10);
}

int main() {
    std::printf("=== Poisson kernel ===\n");
    test_kernel_at_origin();
    test_kernel_known_value();
    test_kernel_symmetry();
    test_kernel_positive_inside();
    test_kernel_boundary_singularity();

    std::printf("=== Poisson integral ===\n");
    test_integral_constant_boundary();
    test_integral_at_origin_equals_mean();
    test_integral_theta_independence_at_origin();

    std::printf("\n%s (%d failure(s))\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
