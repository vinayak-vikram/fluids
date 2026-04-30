// minimal SDL2 viewer for the flapping-airfoil pipeline + freestream + circulation
#include "joukowsky.h"
#include "poisson.h"
#include "velocity_calculator.h"
#include "boundary_builder.h"
#include "lift_calculator.h"

#include <SDL2/SDL.h>
#include <complex>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <cstdio>

using cd = std::complex<double>;

static const int W = 900, H = 650;
static const int GRAPH_H = 140;
static const int SCENE_H = H - GRAPH_H;
static const double VX_MIN = -3.0, VX_MAX = 3.0;
static const double VY_MIN = -2.0, VY_MAX = 2.0;
static const int GX = 90, GY = 60;
static const int N = 48;

static inline void w2px(double wx, double wy, int& px, int& py) {
    px = int((wx - VX_MIN) / (VX_MAX - VX_MIN) * W);
    py = int((VY_MAX - wy) / (VY_MAX - VY_MIN) * SCENE_H);
}

static SDL_Color cmap(double s) {
    s = std::max(-1.0, std::min(1.0, s));
    Uint8 r, g, b;
    if (s >= 0) { r = 255; g = Uint8(255 * (1 - s)); b = g; }
    else        { b = 255; g = Uint8(255 * (1 + s)); r = g; }
    return {r, g, b, 255};
}

int main() {
    const double a = 1.0;
    const cd b{-0.1, 0.08};
    JoukowskyTransform jt(a, b);
    const double r0 = jt.circleRadius();

    cd hinge = jt.z2w(b - r0);
    for (int i = 0; i < 1000; ++i) {
        double th = 2 * M_PI * i / 1000.0;
        cd w = jt.z2w(b + r0 * std::exp(cd(0, th)));
        if (w.real() < hinge.real()) hinge = w;
    }

    PoissonSolver        ps(r0);
    VelocityCalculator   vc(ps, jt);
    BoundaryValueBuilder bb(jt, hinge, N);

    // tunable params
    double U_inf = 2.0;       // freestream in +x (w-plane)
    double Gamma = 0.0;       // circulation, ccw positive
    double flapAmp = 3.0;
    double flapFreq = 0.6;
    bool   flapping = true;
    bool   autoKutta = false;
    const double rho = 1.225;

    // Kutta probe point: just outside trailing edge in z-plane (z=a is the
    // critical point where f'(z)=0). Linear in Gamma → closed-form least-squares.
    cd ea_te = (cd(a, 0.0) - b) / r0;       // unit outward radial at z=a
    cd z_kutta = cd(a, 0.0) + 0.02 * ea_te; // small offset outside wing
    cd B_kutta;
    {
        cd zeta = z_kutta - b;
        cd Q_z_unit = cd(0, -1) / (2.0 * M_PI * zeta);
        B_kutta = Q_z_unit / jt.dz2w(z_kutta);
    }

    // Total complex velocity (Q_w = u-iv) at z (outside airfoil)
    auto Q_w_total = [&](cd z, const std::vector<double>& psi) -> cd {
        cd zeta = z - b;
        cd Q_z_extra = U_inf
                     - U_inf * (r0 * r0) / (zeta * zeta)
                     + cd(0, -1) * Gamma / (2.0 * M_PI * zeta);
        cd Q_w_motion = vc.velocity(z, psi); // already Q_w
        return Q_w_motion + Q_z_extra / jt.dz2w(z);
    };

    // Recompute lift via Bernoulli on total velocity
    auto computeLift = [&](const std::vector<double>& psi) -> double {
        double dtheta = 2 * M_PI / N;
        double rho_eval = r0 * 1.02;
        cd F = 0;
        for (int i = 0; i < N; ++i) {
            double th = i * dtheta;
            cd e_ith = std::exp(cd(0, th));
            cd z_eval = b + rho_eval * e_ith;
            cd q = Q_w_total(z_eval, psi);
            double V2 = std::norm(q);
            cd z_bnd = b + r0 * e_ith;
            cd n_ds  = jt.dz2w(z_bnd) * e_ith * r0 * dtheta;
            F += 0.5 * rho * V2 * n_ds;
        }
        return F.imag();
    };

    // airfoil outline
    std::vector<SDL_Point> airfoil(200);
    for (int i = 0; i < 200; ++i) {
        double th = 2 * M_PI * i / 200.0;
        cd w = jt.z2w(b + r0 * std::exp(cd(0, th)));
        int px, py; w2px(w.real(), w.imag(), px, py);
        airfoil[i] = {px, py};
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window*   win = SDL_CreateWindow("flapping airfoil",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, 0);
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    struct P { double x, y; };
    std::vector<P> parts(180);
    auto resetP = [&](P& p, bool spawnLeft) {
        p.x = spawnLeft ? VX_MIN + 0.05
                        : VX_MIN + (rand() / double(RAND_MAX)) * (VX_MAX - VX_MIN);
        p.y = VY_MIN + (rand() / double(RAND_MAX)) * (VY_MAX - VY_MIN);
    };
    for (auto& p : parts) resetP(p, false);

    std::deque<double> liftHist;
    const int HIST = 400;

    auto setTitle = [&]() {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "airfoil  U=%.2f  Gamma=%.2f%s  flap=%s  [arrows: U/Gamma  K: auto-Kutta  Space: flap  R: reset]",
            U_inf, Gamma, autoKutta ? " [auto]" : "", flapping ? "on" : "off");
        SDL_SetWindowTitle(win, buf);
    };
    setTitle();

    bool running = true;
    Uint32 t0 = SDL_GetTicks();
    Uint32 prev = t0;

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            if (e.type == SDL_KEYDOWN) {
                bool changed = true;
                switch (e.key.keysym.sym) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_UP:     U_inf += 0.5; break;
                    case SDLK_DOWN:   U_inf -= 0.5; break;
                    case SDLK_RIGHT:  Gamma += 0.5; break;
                    case SDLK_LEFT:   Gamma -= 0.5; break;
                    case SDLK_SPACE:  flapping = !flapping; break;
                    case SDLK_k:      autoKutta = !autoKutta; break;
                    case SDLK_r:      U_inf = 2.0; Gamma = 0.0; flapping = true; autoKutta = false; break;
                    default: changed = false;
                }
                if (changed) setTitle();
            }
        }

        Uint32 now = SDL_GetTicks();
        double t  = (now - t0) / 1000.0;
        double dt = (now - prev) / 1000.0;
        prev = now;

        double omega = flapping ? flapAmp * std::sin(2 * M_PI * flapFreq * t) : 0.0;

        auto psi = bb.build(omega);

        if (autoKutta) {
            // Q_w_total = A + Gamma*B at z_kutta; minimize |A+ΓB|² over real Γ.
            double saveG = Gamma; Gamma = 0.0;
            cd A = Q_w_total(z_kutta, psi);
            Gamma = saveG;
            Gamma = -(A.real() * B_kutta.real() + A.imag() * B_kutta.imag())
                    / std::norm(B_kutta);
            setTitle();
        }

        double L = computeLift(psi);
        liftHist.push_back(L);
        if ((int)liftHist.size() > HIST) liftHist.pop_front();

        static int dbgFrame = 0;
        if (++dbgFrame % 30 == 0) {
            std::printf("U=%.2f  Gamma=%+.3f%s  omega=%+.2f  L=%+.4f  KJ_predict=%+.4f\n",
                U_inf, Gamma, autoKutta ? "[K]" : "   ", omega, L, -rho * U_inf * Gamma);
            std::fflush(stdout);
        }

        // ── pressure heatmap ────────────────────────────────────────────
        std::vector<double> V2(GX * GY, 0.0);
        std::vector<bool>   ins(GX * GY, false);
        double V2max = 1e-6;
        for (int gy = 0; gy < GY; ++gy) {
            for (int gx = 0; gx < GX; ++gx) {
                double wx = VX_MIN + (VX_MAX - VX_MIN) * (gx + 0.5) / GX;
                double wy = VY_MAX - (VY_MAX - VY_MIN) * (gy + 0.5) / GY;
                cd w(wx, wy);
                cd z = jt.w2z(w);
                if (std::abs(z - b) < r0 * 1.02) { ins[gy * GX + gx] = true; continue; }
                cd q = Q_w_total(z, psi);
                double v2 = std::norm(q);
                V2[gy * GX + gx] = v2;
                if (v2 < 100.0) V2max = std::max(V2max, v2);
            }
        }

        SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
        SDL_RenderClear(ren);
        double cellW = W / double(GX);
        double cellH = SCENE_H / double(GY);
        for (int gy = 0; gy < GY; ++gy) {
            for (int gx = 0; gx < GX; ++gx) {
                SDL_Rect r{int(gx * cellW), int(gy * cellH),
                           int(cellW + 1), int(cellH + 1)};
                if (ins[gy * GX + gx]) {
                    SDL_SetRenderDrawColor(ren, 40, 40, 40, 255);
                } else {
                    double v2 = std::min(V2[gy * GX + gx], V2max);
                    double s = 1.0 - 2.0 * (v2 / V2max);
                    SDL_Color c = cmap(s);
                    SDL_SetRenderDrawColor(ren, c.r, c.g, c.b, 255);
                }
                SDL_RenderFillRect(ren, &r);
            }
        }

        // ── particles advected by total velocity ─────────────────────────
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 220);
        for (auto& p : parts) {
            cd w(p.x, p.y);
            cd z = jt.w2z(w);
            double vx = U_inf, vy = 0.0;
            if (std::abs(z - b) > r0 * 1.05) {
                cd Q = Q_w_total(z, psi); // u - iv
                cd q = std::conj(Q);
                if (std::isfinite(q.real()) && std::isfinite(q.imag()) && std::abs(q) < 30.0) {
                    vx = q.real(); vy = q.imag();
                }
            }
            p.x += vx * dt;
            p.y += vy * dt;
            if (p.x > VX_MAX || p.x < VX_MIN || p.y > VY_MAX || p.y < VY_MIN
                || std::abs(jt.w2z(cd(p.x, p.y)) - b) < r0 * 1.05) {
                resetP(p, U_inf >= 0);
                if (U_inf < 0) p.x = VX_MAX - 0.05;
            }
            int px, py; w2px(p.x, p.y, px, py);
            if (py >= 0 && py < SCENE_H) {
                SDL_Rect dot{px, py, 2, 2};
                SDL_RenderFillRect(ren, &dot);
            }
        }

        // ── airfoil outline ─────────────────────────────────────────────
        SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
        SDL_RenderDrawLines(ren, airfoil.data(), (int)airfoil.size());
        SDL_RenderDrawLine(ren, airfoil.back().x, airfoil.back().y,
                           airfoil.front().x, airfoil.front().y);

        int hx, hy; w2px(hinge.real(), hinge.imag(), hx, hy);
        SDL_SetRenderDrawColor(ren, 255, 255, 0, 255);
        SDL_Rect hr{hx - 3, hy - 3, 6, 6};
        SDL_RenderFillRect(ren, &hr);

        // ── lift graph ──────────────────────────────────────────────────
        SDL_Rect graph{0, SCENE_H, W, GRAPH_H};
        SDL_SetRenderDrawColor(ren, 20, 20, 25, 255);
        SDL_RenderFillRect(ren, &graph);

        double Lmax = 1.0;
        for (double l : liftHist) Lmax = std::max(Lmax, std::abs(l));
        int midY = SCENE_H + GRAPH_H / 2;
        SDL_SetRenderDrawColor(ren, 80, 80, 80, 255);
        SDL_RenderDrawLine(ren, 0, midY, W, midY);

        int prevX = -1, prevY = 0;
        for (size_t i = 0; i < liftHist.size(); ++i) {
            int gxp = int(i * W / double(HIST));
            int gyp = midY - int(liftHist[i] / Lmax * (GRAPH_H / 2 - 8));
            if (prevX >= 0) {
                if (liftHist[i] >= 0) SDL_SetRenderDrawColor(ren, 100, 255, 100, 255);
                else                  SDL_SetRenderDrawColor(ren, 255, 80, 80, 255);
                SDL_RenderDrawLine(ren, prevX, prevY, gxp, gyp);
            }
            prevX = gxp; prevY = gyp;
        }

        int curY = midY - int(L / Lmax * (GRAPH_H / 2 - 8));
        if (L >= 0) SDL_SetRenderDrawColor(ren, 100, 255, 100, 255);
        else        SDL_SetRenderDrawColor(ren, 255, 80, 80, 255);
        SDL_Rect bar{W - 6, std::min(curY, midY), 6, std::abs(curY - midY) + 1};
        SDL_RenderFillRect(ren, &bar);

        SDL_RenderPresent(ren);
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
