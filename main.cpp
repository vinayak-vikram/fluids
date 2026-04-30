// minimal SDL2 viewer: stationary Joukowsky airfoil in freestream U with circulation Gamma
#include "joukowsky.h"

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
static const int GX = 110, GY = 75;

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
    const double rho = 1.225;

    double U_inf = 2.0;
    double Gamma = 0.0;
    bool   autoKutta = false;

    // Q_w_total at z (just outside airfoil): freestream + dipole + free vortex.
    auto Q_w_total = [&](cd z) -> cd {
        cd zeta = z - b;
        cd Q_z = U_inf
               - U_inf * (r0 * r0) / (zeta * zeta)
               + cd(0, -1) * Gamma / (2.0 * M_PI * zeta);
        return Q_z / jt.dz2w(z);
    };

    const int N = 400;
    auto computeLift = [&]() -> double {
        double dtheta = 2 * M_PI / N;
        double rho_eval = r0 * 1.02;
        cd F = 0;
        for (int i = 0; i < N; ++i) {
            double th = i * dtheta;
            cd e_ith = std::exp(cd(0, th));
            cd z_eval = b + rho_eval * e_ith;
            double V2 = std::norm(Q_w_total(z_eval));
            cd z_bnd = b + r0 * e_ith;
            cd n_ds  = jt.dz2w(z_bnd) * e_ith * r0 * dtheta;
            F += 0.5 * rho * V2 * n_ds;
        }
        return F.imag();
    };

    // Kutta probe: just outside z=a in radial direction; Q_w linear in Gamma.
    cd ea_te   = (cd(a, 0.0) - b) / r0;
    cd z_kutta = cd(a, 0.0) + 0.02 * ea_te;
    cd B_kutta;
    {
        cd zeta = z_kutta - b;
        B_kutta = (cd(0, -1) / (2.0 * M_PI * zeta)) / jt.dz2w(z_kutta);
    }
    auto applyKutta = [&]() {
        double saveG = Gamma; Gamma = 0.0;
        cd A = Q_w_total(z_kutta);
        Gamma = saveG;
        Gamma = -(A.real() * B_kutta.real() + A.imag() * B_kutta.imag())
                / std::norm(B_kutta);
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
    SDL_Window*   win = SDL_CreateWindow("airfoil",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W, H, 0);
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    struct P { double x, y; };
    std::vector<P> parts(220);
    auto resetP = [&](P& p) {
        p.x = U_inf >= 0 ? VX_MIN + 0.05 : VX_MAX - 0.05;
        p.y = VY_MIN + (rand() / double(RAND_MAX)) * (VY_MAX - VY_MIN);
    };
    for (auto& p : parts) {
        p.x = VX_MIN + (rand() / double(RAND_MAX)) * (VX_MAX - VX_MIN);
        p.y = VY_MIN + (rand() / double(RAND_MAX)) * (VY_MAX - VY_MIN);
    }

    std::deque<double> liftHist;
    const int HIST = 400;

    auto setTitle = [&]() {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "airfoil  U=%.2f  Gamma=%+.3f%s  L=%+.3f  [arrows: U/Gamma  K: auto-Kutta  R: reset]",
            U_inf, Gamma, autoKutta ? " [auto]" : "",
            liftHist.empty() ? 0.0 : liftHist.back());
        SDL_SetWindowTitle(win, buf);
    };
    setTitle();

    bool running = true;
    Uint32 prev = SDL_GetTicks();

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_UP:     U_inf += 0.5; break;
                    case SDLK_DOWN:   U_inf -= 0.5; break;
                    case SDLK_RIGHT:  Gamma += 0.5; break;
                    case SDLK_LEFT:   Gamma -= 0.5; break;
                    case SDLK_k:      autoKutta = !autoKutta; break;
                    case SDLK_r:      U_inf = 2.0; Gamma = 0.0; autoKutta = false; break;
                    default: break;
                }
            }
        }

        Uint32 now = SDL_GetTicks();
        double dt = (now - prev) / 1000.0;
        prev = now;

        if (autoKutta) applyKutta();

        double L = computeLift();
        liftHist.push_back(L);
        if ((int)liftHist.size() > HIST) liftHist.pop_front();
        setTitle();

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
                double v2 = std::norm(Q_w_total(z));
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
                cd q = std::conj(Q_w_total(z));
                if (std::isfinite(q.real()) && std::isfinite(q.imag()) && std::abs(q) < 30.0) {
                    vx = q.real(); vy = q.imag();
                }
            }
            p.x += vx * dt;
            p.y += vy * dt;
            if (p.x > VX_MAX || p.x < VX_MIN || p.y > VY_MAX || p.y < VY_MIN
                || std::abs(jt.w2z(cd(p.x, p.y)) - b) < r0 * 1.05) {
                resetP(p);
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
