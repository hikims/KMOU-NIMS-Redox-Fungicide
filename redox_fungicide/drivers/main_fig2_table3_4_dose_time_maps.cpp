/*
 * Driver for the main-manuscript Table 4 threshold comparison.
 * Historical filename retained.
 */
#include "functions.hpp"
#include "initial.hpp"
#include "parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <chrono>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#ifndef OUTPUT_ROOT
#define OUTPUT_ROOT "./results/"
#endif
#ifndef PATH
#define PATH OUTPUT_ROOT
#endif

#ifndef MACHINE_LABEL
#define MACHINE_LABEL "macstudio"
#endif

// ============================================================================
// MATCOM Tier-B: B2-only numerical driver
//
// Purpose:
//   Quantitative comparison of analytical sufficient thresholds and empirical
//   fixed-time thresholds as alpha varies.
//
//   Revised after B1: the fixed observation time is Tobs=12 and the
//   empirical criteria are Mbar < 2.5e-2 and M_inf < 2.0e-1.
//
// Output directory:
//   PATH/TierB/summary/
//
// Required case file:
//   PATH/TierB/b2_cases.txt
// with columns:
//   case_id alpha_factor A0_min A0_max
// ============================================================================

// ----------------------------------------------------------------------------
// USER CONFIGURATION
// ----------------------------------------------------------------------------
static const double MSTAR = 2.5e-2;   // threshold for normalized aggregate burden Mbar at Tobs=12
static const double FSTAR = 2.0e-1;   // threshold for L-infinity burden at Tobs=12

static const double B2_TOBS          = 12.0;
static const double B2_BETA          = 8.0;   // persistence regime for theta=0.3, m=2
static const int    B2_BISECT_ITERS  = 8;     // bisection iterations after endpoint bracketing

// Analytical threshold: numerical minimization over admissible A_aux.
static const int    Aaux_GRID_N      = 5000;
static const double Aaux_MIN_FACTOR  = 1.0001;
static const double Aaux_MAX         = 500.0;

// Numerical output control.
static const bool   VERBOSE_EACH_RUN = true;

// ----------------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------------
static void mkdir_p(const std::string& path) {
    if (path.empty()) return;
    std::string cur;
    size_t i = 0;
    if (path[0] == '/') {
        cur = "/";
        i = 1;
    }
    while (i <= path.size()) {
        size_t j = path.find('/', i);
        std::string part = path.substr(i, j == std::string::npos ? std::string::npos : j - i);
        if (!part.empty()) {
            if (!cur.empty() && cur.back() != '/') cur += "/";
            cur += part;
            struct stat st;
            if (stat(cur.c_str(), &st) == -1) {
                if (mkdir(cur.c_str(), 0755) != 0 && errno != EEXIST) {
                    perror(("mkdir failed: " + cur).c_str());
                    exit(1);
                }
            } else if (!S_ISDIR(st.st_mode)) {
                fprintf(stderr, "Path exists but is not a directory: %s\n", cur.c_str());
                exit(1);
            }
        }
        if (j == std::string::npos) break;
        i = j + 1;
    }
}

static void make_output_dirs() {
    mkdir_p(std::string(PATH));
    mkdir_p(std::string(PATH) + "TierB");
    mkdir_p(std::string(PATH) + "TierB/summary");
    mkdir_p(std::string(PATH) + "TierB/logs");
}

static double safe_log10(double x) {
    if (x <= 0.0 || !std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
    return std::log10(x);
}

static double domain_area() {
    return Lx * Ly;
}

static double fungal_area_percent(double **F) {
    int infected_cells = 0;
    int total_cells = (Nx - 1) * (Ny - 1);

    for (int i = 0; i < Nx - 1; i++) {
        for (int j = 0; j < Ny - 1; j++) {
            double Fav = 0.25 * (F[i][j] + F[i+1][j] + F[i][j+1] + F[i+1][j+1]);
            if (Fav >= F_area_threshold) infected_cells++;
        }
    }
    return 100.0 * (double)infected_cells / (double)total_cells;
}

static double compute_total_mass(double **U) {
    int a = 0, b = Nx - 1, c = 0, d = Ny - 1;
    double first_term = 0.0, second_term = 0.0;
    double third_term = 0.0, fourth_term = 0.0, fifth_term = 0.0;

    for (int j = 1; j < Ny - 1; j++) {
        first_term += U[a][j];
        second_term += U[b][j];
    }
    for (int i = 1; i < Nx - 1; i++) {
        third_term += U[i][c];
        fourth_term += U[i][d];
    }
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            fifth_term += U[i][j];
        }
    }

    return 0.25 * dx * dy *
           (U[a][c] + U[a][d] + U[b][c] + U[b][d]
            + 2.0 * (first_term + second_term + third_term + fourth_term)
            + 4.0 * fifth_term);
}

static void minmax_field(double **U, double &umin, double &umax) {
    umin =  std::numeric_limits<double>::infinity();
    umax = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            umin = std::min(umin, U[i][j]);
            umax = std::max(umax, U[i][j]);
        }
    }
}

// Local Hill functions with runtime beta.
static inline double hill_f_run(double Rval, double m_run, double theta_run) {
    double Rm = std::pow(Rval, m_run);
    double thm = std::pow(theta_run, m_run);
    return Rm / (Rm + thm);
}

static inline double hill_g_run(double Rval, double beta_run, double m_run, double theta_run) {
    double Rm = std::pow(Rval, m_run);
    double thm = std::pow(theta_run, m_run);
    return beta_run * thm / (Rm + thm);
}

static inline double coefficient_dF_beta(double Fij,
                                         double r,
                                         double Fp,
                                         double Fc,
                                         double Fm,
                                         double Rval,
                                         double m_run,
                                         double theta_run,
                                         double beta_run) {
    double fR = hill_f_run(Rval, m_run, theta_run);
    double gR = hill_g_run(Rval, beta_run, m_run, theta_run);
    double reaction = fR * Fij * (1.0 - Fij) - gR * Fij;
    return Fij + r * (Fp - 2.0 * Fc + Fm) + dt * 0.5 * reaction;
}

struct Observation {
    double t;
    double total_A;
    double total_F;
    double Mbar;
    double M_inf;
    double area_percent;
    double min_A, max_A;
    double min_R, max_R;
    double min_F, max_F;
};

struct RunResult {
    std::string tag;
    double A0;
    double beta_run;
    double dA_run;
    double dF_run;
    double alpha_run;
    Observation final_obs;
    bool unstable;
    double wall_seconds;
};

struct CaseDef {
    std::string case_id;
    double alpha_factor;
    double A0_min;
    double A0_max;
};

struct EvalRecord {
    std::string case_id;
    double alpha_factor;
    double alpha_run;
    double A0;
    RunResult result;
};

static Observation collect_observation(double t, double **A, double **F, double **R) {
    Observation o;
    o.t = t;
    o.total_A = compute_total_mass(A);
    o.total_F = compute_total_mass(F);
    o.Mbar = o.total_F / domain_area();
    o.area_percent = fungal_area_percent(F);

    minmax_field(A, o.min_A, o.max_A);
    minmax_field(R, o.min_R, o.max_R);
    minmax_field(F, o.min_F, o.max_F);
    o.M_inf = o.max_F;
    return o;
}

static RunResult run_simulation(
    const std::string& tag,
    double A0,
    double beta_run,
    double dA_run,
    double dF_run,
    double alpha_run,
    double t_final_requested
) {
    auto wall0 = std::chrono::steady_clock::now();

    RunResult result;
    result.tag = tag;
    result.A0 = A0;
    result.beta_run = beta_run;
    result.dA_run = dA_run;
    result.dF_run = dF_run;
    result.alpha_run = alpha_run;
    result.unstable = false;
    result.wall_seconds = 0.0;

    double **A = mallocDouble(Nx, Ny);
    double **R = mallocDouble(Nx, Ny);
    double **F = mallocDouble(Nx, Ny);
    double **A_new = mallocDouble(Nx, Ny);
    double **R_new = mallocDouble(Nx, Ny);
    double **F_new = mallocDouble(Nx, Ny);

    double *A_Ax = mallocSingle(Nx);
    double *B_Ax = mallocSingle(Nx);
    double *C_Ax = mallocSingle(Nx);
    double *A_Ay = mallocSingle(Ny);
    double *B_Ay = mallocSingle(Ny);
    double *C_Ay = mallocSingle(Ny);
    double **x_Ax = mallocDouble(Ny, Nx);
    double **x_Ay = mallocDouble(Nx, Ny);
    double **z_Ax = mallocDouble(Ny, Nx);
    double **z_Ay = mallocDouble(Nx, Ny);

    double *A_Fx = mallocSingle(Nx);
    double *B_Fx = mallocSingle(Nx);
    double *C_Fx = mallocSingle(Nx);
    double *A_Fy = mallocSingle(Ny);
    double *B_Fy = mallocSingle(Ny);
    double *C_Fy = mallocSingle(Ny);
    double **x_Fx = mallocDouble(Ny, Nx);
    double **x_Fy = mallocDouble(Nx, Ny);
    double **z_Fx = mallocDouble(Ny, Nx);
    double **z_Fy = mallocDouble(Nx, Ny);

    double **Ax = mallocDouble(Ny, Nx);
    double **Ay = mallocDouble(Nx, Ny);
    double **Fx = mallocDouble(Ny, Nx);
    double **Fy = mallocDouble(Nx, Ny);

    double **L_Ax = mallocDouble(Lp, Nx);
    double **L_Ay = mallocDouble(Lp, Ny);
    double *U_Ax = mallocSingle(Nx);
    double *U_Ay = mallocSingle(Ny);
    double **L_Fx = mallocDouble(Lp, Nx);
    double **L_Fy = mallocDouble(Lp, Ny);
    double *U_Fx = mallocSingle(Nx);
    double *U_Fy = mallocSingle(Ny);

    double rAx = (dA_run * 0.5 * dt) / (dx * dx);
    double rAy = (dA_run * 0.5 * dt) / (dy * dy);
    double rFx = (dF_run * 0.5 * dt) / (dx * dx);
    double rFy = (dF_run * 0.5 * dt) / (dy * dy);

    srand(1);
    initialization(R, Nx, Ny, initR);
    initialization(A, Nx, Ny, A0);
    Circular(F, Nx, Ny, dx, dy,
             F_center_cx, F_center_cy,
             F_radius, F_inside, F_outside);

    for (int i = 0; i < Nx; i++) { A_Ax[i] = 1 + 2*rAx; B_Ax[i] = -rAx; C_Ax[i] = -rAx; }
    A_Ax[0] = 1 + 2*rAx; B_Ax[0] = 0.0; C_Ax[0] = -2*rAx;
    A_Ax[Nx-1] = 1 + 2*rAx; B_Ax[Nx-1] = -2*rAx; C_Ax[Nx-1] = 0.0;

    for (int j = 0; j < Ny; j++) { A_Ay[j] = 1 + 2*rAy; B_Ay[j] = -rAy; C_Ay[j] = -rAy; }
    A_Ay[0] = 1 + 2*rAy; B_Ay[0] = 0.0; C_Ay[0] = -2*rAy;
    A_Ay[Ny-1] = 1 + 2*rAy; B_Ay[Ny-1] = -2*rAy; C_Ay[Ny-1] = 0.0;
    LUx(A_Ax, B_Ax, C_Ax, L_Ax, U_Ax, Nx);
    LUy(A_Ay, B_Ay, C_Ay, L_Ay, U_Ay, Ny);

    for (int i = 0; i < Nx; i++) { A_Fx[i] = 1 + 2*rFx; B_Fx[i] = -rFx; C_Fx[i] = -rFx; }
    A_Fx[0] = 1 + 2*rFx; B_Fx[0] = 0.0; C_Fx[0] = -2*rFx;
    A_Fx[Nx-1] = 1 + 2*rFx; B_Fx[Nx-1] = -2*rFx; C_Fx[Nx-1] = 0.0;

    for (int j = 0; j < Ny; j++) { A_Fy[j] = 1 + 2*rFy; B_Fy[j] = -rFy; C_Fy[j] = -rFy; }
    A_Fy[0] = 1 + 2*rFy; B_Fy[0] = 0.0; C_Fy[0] = -2*rFy;
    A_Fy[Ny-1] = 1 + 2*rFy; B_Fy[Ny-1] = -2*rFy; C_Fy[Ny-1] = 0.0;
    LUx(A_Fx, B_Fx, C_Fx, L_Fx, U_Fx, Nx);
    LUy(A_Fy, B_Fy, C_Fy, L_Fy, U_Fy, Ny);

    int max_steps_requested = (int)std::ceil(t_final_requested / dt);
    int max_steps = std::min(iter, max_steps_requested);
    if (max_steps < max_steps_requested) {
        printf("Warning [%s]: requested final time %.6lf exceeds iter*dt %.6lf. Running to %.6lf only.\n",
               tag.c_str(), t_final_requested, iter * dt, max_steps * dt);
    }

    tbb::task_arena arena;

    for (int k = 0; k < max_steps; k++) {
        for (int j = 0; j < Ny; j++) {
            A[0][j]      = A[1][j];
            A[Nx-1][j]   = A[Nx-2][j];
            F[0][j]      = F[1][j];
            F[Nx-1][j]   = F[Nx-2][j];
        }
        for (int i = 0; i < Nx; i++) {
            A[i][0]      = A[i][1];
            A[i][Ny-1]   = A[i][Ny-2];
            F[i][0]      = F[i][1];
            F[i][Ny-1]   = F[i][Ny-2];
        }

        arena.execute([&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, Ny),
                [&](const tbb::blocked_range<int> &r) {
                    for (int j = r.begin(); j < r.end(); j++) {
                        for (int i = 0; i < Nx; i++) {
                            int jp = (j == Ny - 1) ? j - 1 : j + 1;
                            int jm = (j == 0)      ? j + 1 : j - 1;
                            Ax[j][i] = coefficient_dA(A[i][j], rAy,
                                                      A[i][jp], A[i][j], A[i][jm],
                                                      F[i][j], alpha_run);
                            Fx[j][i] = coefficient_dF_beta(F[i][j], rFy,
                                                           F[i][jp], F[i][j], F[i][jm],
                                                           R[i][j], m, theta, beta_run);
                        }
                        solx(L_Ax, U_Ax, &Ax[j][0], &z_Ax[j][0], &x_Ax[j][0], Nx);
                        solx(L_Fx, U_Fx, &Fx[j][0], &z_Fx[j][0], &x_Fx[j][0], Nx);
                        for (int i = 0; i < Nx; i++) {
                            A_new[i][j] = x_Ax[j][i];
                            F_new[i][j] = x_Fx[j][i];
                            R_new[i][j] = R[i][j] + dt * (eps * (1.0 - R[i][j]) - A[i][j] * R[i][j]);
                        }
                    }
                });
        });

        arena.execute([&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, Nx),
                [&](const tbb::blocked_range<int> &r) {
                    for (int i = r.begin(); i < r.end(); i++) {
                        for (int j = 0; j < Ny; j++) {
                            int ip = (i == Nx - 1) ? i - 1 : i + 1;
                            int im = (i == 0)      ? i + 1 : i - 1;
                            Ay[i][j] = coefficient_dA(A_new[i][j], rAx,
                                                      A_new[ip][j], A_new[i][j], A_new[im][j],
                                                      F_new[i][j], alpha_run);
                            Fy[i][j] = coefficient_dF_beta(F_new[i][j], rFx,
                                                           F_new[ip][j], F_new[i][j], F_new[im][j],
                                                           R_new[i][j], m, theta, beta_run);
                        }
                        soly(L_Ay, U_Ay, &Ay[i][0], &z_Ay[i][0], &x_Ay[i][0], Ny);
                        soly(L_Fy, U_Fy, &Fy[i][0], &z_Fy[i][0], &x_Fy[i][0], Ny);
                        for (int j = 0; j < Ny; j++) {
                            A[i][j] = x_Ay[i][j];
                            if (A[i][j] < 0.0) A[i][j] = 0.0;
                            F[i][j] = x_Fy[i][j];
                            if (F[i][j] < 0.0) F[i][j] = 0.0;
                            if (F[i][j] > 1.0) F[i][j] = 1.0;
                            R[i][j] = R_new[i][j];
                            if (R[i][j] < 0.0) R[i][j] = 0.0;
                            if (R[i][j] > 1.0) R[i][j] = 1.0;
                        }
                    }
                });
        });

        if (!std::isfinite(F[Nx/2][Ny/2]) || !std::isfinite(A[Nx/2][Ny/2]) || !std::isfinite(R[Nx/2][Ny/2])) {
            printf("Unstable at step %d (t = %.6lf).\n", k+1, (k+1)*dt);
            result.unstable = true;
            break;
        }
    }

    result.final_obs = collect_observation(max_steps * dt, A, F, R);

    freeMallocDouble(L_Ax, Lp); freeMallocDouble(L_Ay, Lp);
    freeMallocDouble(L_Fx, Lp); freeMallocDouble(L_Fy, Lp);
    freeMallocSingle(U_Ax); freeMallocSingle(U_Ay);
    freeMallocSingle(U_Fx); freeMallocSingle(U_Fy);

    freeMallocDouble(A, Nx); freeMallocDouble(A_new, Nx);
    freeMallocDouble(R, Nx); freeMallocDouble(R_new, Nx);
    freeMallocDouble(F, Nx); freeMallocDouble(F_new, Nx);

    freeMallocDouble(x_Ax, Ny); freeMallocDouble(x_Ay, Nx);
    freeMallocDouble(z_Ax, Ny); freeMallocDouble(z_Ay, Nx);
    freeMallocSingle(A_Ax); freeMallocSingle(B_Ax); freeMallocSingle(C_Ax);
    freeMallocSingle(A_Ay); freeMallocSingle(B_Ay); freeMallocSingle(C_Ay);

    freeMallocDouble(x_Fx, Ny); freeMallocDouble(x_Fy, Nx);
    freeMallocDouble(z_Fx, Ny); freeMallocDouble(z_Fy, Nx);
    freeMallocSingle(A_Fx); freeMallocSingle(B_Fx); freeMallocSingle(C_Fx);
    freeMallocSingle(A_Fy); freeMallocSingle(B_Fy); freeMallocSingle(C_Fy);

    freeMallocDouble(Ax, Ny); freeMallocDouble(Ay, Nx);
    freeMallocDouble(Fx, Ny); freeMallocDouble(Fy, Nx);

    auto wall1 = std::chrono::steady_clock::now();
    result.wall_seconds = std::chrono::duration<double>(wall1 - wall0).count();

    return result;
}

// ----------------------------------------------------------------------------
// Analytical sufficient threshold from Theorem 4.3, computed in log form.
// ----------------------------------------------------------------------------
static double log_A0_suf_given_Aaux(double Aaux,
                                    double alpha_run,
                                    double beta_run,
                                    double F_star) {
    if (beta_run <= 1.0) return std::numeric_limits<double>::infinity();
    double theta_m = std::pow(theta, m);
    if (beta_run * theta_m >= 1.0) return std::numeric_limits<double>::infinity();

    double Aaux_min = eps * (1.0 - theta) / theta;
    if (Aaux <= Aaux_min) return std::numeric_limits<double>::infinity();

    double r_aux = eps / (eps + Aaux);
    if (r_aux >= theta) return std::numeric_limits<double>::infinity();

    double tin = (1.0 / (eps + Aaux)) * std::log((1.0 - r_aux) / (theta - r_aux));
    if (!std::isfinite(tin) || tin < 0.0) return std::numeric_limits<double>::infinity();

    double M0 = std::min(1.0, std::max(F_inside, F_outside));
    double Mpre = std::min(1.0, M0 * std::exp(tin));
    if (F_star <= 0.0 || F_star >= Mpre) return std::numeric_limits<double>::infinity();

    double mu0 = 0.5 * (beta_run - 1.0);
    double lambda = alpha_run + l;

    double log_A0 = std::log(Aaux) + lambda * (tin + (1.0/mu0) * std::log(Mpre / F_star));
    return log_A0;
}

static void optimize_A0_suf(double alpha_run,
                            double beta_run,
                            double F_star,
                            double &best_Aaux,
                            double &best_logA0,
                            double &best_log10A0) {
    best_Aaux = std::numeric_limits<double>::quiet_NaN();
    best_logA0 = std::numeric_limits<double>::infinity();
    best_log10A0 = std::numeric_limits<double>::infinity();

    double Aaux_min = eps * (1.0 - theta) / theta * Aaux_MIN_FACTOR;
    double Aaux_max = Aaux_MAX;
    if (Aaux_max <= Aaux_min) Aaux_max = Aaux_min * 100.0;

    double lmin = std::log(Aaux_min);
    double lmax = std::log(Aaux_max);
    for (int i = 0; i < Aaux_GRID_N; i++) {
        double s = (Aaux_GRID_N == 1) ? 0.0 : (double)i / (double)(Aaux_GRID_N - 1);
        double Aaux = std::exp(lmin + s * (lmax - lmin));
        double logA0 = log_A0_suf_given_Aaux(Aaux, alpha_run, beta_run, F_star);
        if (std::isfinite(logA0) && logA0 < best_logA0) {
            best_logA0 = logA0;
            best_Aaux = Aaux;
        }
    }
    best_log10A0 = best_logA0 / std::log(10.0);
}

static bool ok_Mbar(const RunResult& r) { return r.final_obs.Mbar  < MSTAR; }
static bool ok_Minf(const RunResult& r) { return r.final_obs.M_inf < FSTAR; }

static RunResult evaluate_A0(const CaseDef& c,
                             double A0,
                             std::vector<EvalRecord>& cache,
                             FILE* fScan) {
    for (const auto& e : cache) {
        if (e.case_id == c.case_id && std::fabs(e.A0 - A0) < 1.0e-12) {
            return e.result;
        }
    }

    double alpha_run = alpha * c.alpha_factor;
    char tag[256];
    snprintf(tag, sizeof(tag), "%s_A0_%.8f", c.case_id.c_str(), A0);

    RunResult r = run_simulation(tag, A0, B2_BETA, dA, dF, alpha_run, B2_TOBS);

    EvalRecord rec;
    rec.case_id = c.case_id;
    rec.alpha_factor = c.alpha_factor;
    rec.alpha_run = alpha_run;
    rec.A0 = A0;
    rec.result = r;
    cache.push_back(rec);

    if (fScan) {
        const Observation& o = r.final_obs;
        fprintf(fScan,
            "%s %s %.12lf %.12lf %.12lf %.12lf %.12lf %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12lf %d %d %.12lf\n",
            MACHINE_LABEL,
            c.case_id.c_str(),
            c.alpha_factor,
            alpha_run,
            A0,
            B2_BETA,
            o.t,
            o.total_A,
            o.total_F,
            o.Mbar,
            o.M_inf,
            o.area_percent,
            o.min_A,
            o.max_A,
            o.min_R,
            o.max_R,
            o.min_F,
            o.max_F,
            safe_log10(o.Mbar),
            (int)ok_Mbar(r),
            (int)ok_Minf(r),
            r.wall_seconds);
        fflush(fScan);
    }

    if (VERBOSE_EACH_RUN) {
        printf("[%s] case=%s A0=%.6lf alpha_factor=%.3lf alpha=%.6lf Mbar=%.6e Minf=%.6e okM=%d okInf=%d wall=%.1fs\n",
               MACHINE_LABEL, c.case_id.c_str(), A0, c.alpha_factor, alpha_run,
               r.final_obs.Mbar, r.final_obs.M_inf, (int)ok_Mbar(r), (int)ok_Minf(r), r.wall_seconds);
        fflush(stdout);
    }
    return r;
}

struct ThresholdResult {
    bool found;
    bool lower_already_ok;
    bool upper_not_ok;
    double lower;
    double upper;
    double mid;
    double value_mid;
};

static ThresholdResult bisection_threshold(const CaseDef& c,
                                           bool use_Minf,
                                           std::vector<EvalRecord>& cache,
                                           FILE* fScan) {
    ThresholdResult tr;
    tr.found = false;
    tr.lower_already_ok = false;
    tr.upper_not_ok = false;
    tr.lower = c.A0_min;
    tr.upper = c.A0_max;
    tr.mid = std::numeric_limits<double>::quiet_NaN();
    tr.value_mid = std::numeric_limits<double>::quiet_NaN();

    RunResult rlo = evaluate_A0(c, c.A0_min, cache, fScan);
    bool lo_ok = use_Minf ? ok_Minf(rlo) : ok_Mbar(rlo);
    if (lo_ok) {
        tr.found = true;
        tr.lower_already_ok = true;
        tr.lower = c.A0_min;
        tr.upper = c.A0_min;
        tr.mid = c.A0_min;
        tr.value_mid = use_Minf ? rlo.final_obs.M_inf : rlo.final_obs.Mbar;
        return tr;
    }

    RunResult rhi = evaluate_A0(c, c.A0_max, cache, fScan);
    bool hi_ok = use_Minf ? ok_Minf(rhi) : ok_Mbar(rhi);
    if (!hi_ok) {
        tr.found = false;
        tr.upper_not_ok = true;
        tr.lower = c.A0_min;
        tr.upper = c.A0_max;
        tr.mid = std::numeric_limits<double>::quiet_NaN();
        tr.value_mid = use_Minf ? rhi.final_obs.M_inf : rhi.final_obs.Mbar;
        return tr;
    }

    double lo = c.A0_min;
    double hi = c.A0_max;
    for (int it = 0; it < B2_BISECT_ITERS; it++) {
        double mid = 0.5 * (lo + hi);
        RunResult rm = evaluate_A0(c, mid, cache, fScan);
        bool ok = use_Minf ? ok_Minf(rm) : ok_Mbar(rm);
        if (ok) hi = mid;
        else    lo = mid;
    }

    tr.found = true;
    tr.lower = lo;
    tr.upper = hi;
    tr.mid = 0.5 * (lo + hi);

    RunResult rmid = evaluate_A0(c, tr.mid, cache, fScan);
    tr.value_mid = use_Minf ? rmid.final_obs.M_inf : rmid.final_obs.Mbar;
    return tr;
}

static std::vector<CaseDef> read_cases() {
    std::vector<CaseDef> cases;
    std::string filename = std::string(PATH) + "TierB/b2_cases.txt";
    std::ifstream fin("configs/b2_cases.txt");
    if (!fin.is_open()) fin.open(filename.c_str());
    if (!fin.is_open()) {
        // Fallback to local working directory.
        fin.open("b2_cases.txt");
    }
    if (!fin.is_open()) {
        fprintf(stderr, "Cannot open b2_cases.txt. Expected at %s or current directory.\n", filename.c_str());
        exit(1);
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        std::istringstream iss(line);
        CaseDef c;
        if (!(iss >> c.case_id >> c.alpha_factor >> c.A0_min >> c.A0_max)) {
            std::cerr << "Skipping invalid case line: " << line << std::endl;
            continue;
        }
        cases.push_back(c);
    }
    return cases;
}

static void run_B2_only() {
    std::vector<CaseDef> cases = read_cases();
    if (cases.empty()) {
        fprintf(stderr, "No B2 cases found.\n");
        exit(1);
    }

    std::string scan_file = std::string(PATH) + "TierB/summary/B2_scan_" + MACHINE_LABEL + ".txt";
    std::string cmp_file  = std::string(PATH) + "TierB/summary/B2_threshold_comparison_" + MACHINE_LABEL + ".txt";

    FILE* fScan = fopen(scan_file.c_str(), "w");
    if (!fScan) { perror("failed to open B2 scan file"); exit(1); }
    FILE* fCmp = fopen(cmp_file.c_str(), "w");
    if (!fCmp) { perror("failed to open B2 comparison file"); fclose(fScan); exit(1); }

    fprintf(fScan,
        "# machine case_id alpha_factor alpha_run A0 beta t total_A total_F Mbar M_inf area_percent min_A max_A min_R max_R min_F max_F log10_Mbar ok_Mbar ok_Minf wall_seconds\n");

    fprintf(fCmp,
        "# machine case_id alpha_factor alpha_run A0_min A0_max "
        "found_Mbar lower_Mbar upper_Mbar mid_Mbar value_Mbar_mid "
        "found_Minf lower_Minf upper_Minf mid_Minf value_Minf_mid "
        "Aaux_opt log10_A0_suf_opt log10_conservativeness_vs_Minf log10_conservativeness_vs_Mbar "
        "endpoint_note_Mbar endpoint_note_Minf n_unique_runs\n");

    for (const auto& c : cases) {
        printf("\n=====================================\n");
        printf("B2 case %s: alpha_factor=%.6lf, A0 range [%.6lf, %.6lf]\n",
               c.case_id.c_str(), c.alpha_factor, c.A0_min, c.A0_max);
        printf("=====================================\n");
        fflush(stdout);

        std::vector<EvalRecord> cache;
        double alpha_run = alpha * c.alpha_factor;

        ThresholdResult trM = bisection_threshold(c, false, cache, fScan);
        ThresholdResult trI = bisection_threshold(c, true,  cache, fScan);

        double Aaux_opt, logA0suf, log10A0suf;
        optimize_A0_suf(alpha_run, B2_BETA, FSTAR, Aaux_opt, logA0suf, log10A0suf);

        double log10_cons_Minf = std::numeric_limits<double>::quiet_NaN();
        double log10_cons_Mbar = std::numeric_limits<double>::quiet_NaN();
        if (trI.found && trI.mid > 0.0 && std::isfinite(log10A0suf)) {
            log10_cons_Minf = log10A0suf - std::log10(trI.mid);
        }
        if (trM.found && trM.mid > 0.0 && std::isfinite(log10A0suf)) {
            log10_cons_Mbar = log10A0suf - std::log10(trM.mid);
        }

        const char* noteM = trM.lower_already_ok ? "lower_already_ok" : (trM.upper_not_ok ? "upper_not_ok" : "bracketed");
        const char* noteI = trI.lower_already_ok ? "lower_already_ok" : (trI.upper_not_ok ? "upper_not_ok" : "bracketed");

        fprintf(fCmp,
            "%s %s %.12lf %.12lf %.12lf %.12lf "
            "%d %.12lf %.12lf %.12lf %.12e "
            "%d %.12lf %.12lf %.12lf %.12e "
            "%.12lf %.12lf %.12lf %.12lf %s %s %zu\n",
            MACHINE_LABEL, c.case_id.c_str(), c.alpha_factor, alpha_run, c.A0_min, c.A0_max,
            (int)trM.found, trM.lower, trM.upper, trM.mid, trM.value_mid,
            (int)trI.found, trI.lower, trI.upper, trI.mid, trI.value_mid,
            Aaux_opt, log10A0suf, log10_cons_Minf, log10_cons_Mbar,
            noteM, noteI, cache.size());
        fflush(fCmp);

        printf("B2 case done: %s | Mbar mid=%.6lf (%s), Minf mid=%.6lf (%s), log10(A0_suf)=%.6lf, unique runs=%zu\n",
               c.case_id.c_str(), trM.mid, noteM, trI.mid, noteI, log10A0suf, cache.size());
        fflush(stdout);
    }

    fclose(fScan);
    fclose(fCmp);

    printf("\nB2 outputs written to:\n%s\n%s\n", scan_file.c_str(), cmp_file.c_str());
}

int main() {
    printf("\n=====================================\n");
    printf("Running MATCOM Tier-B B2-only driver (%s)\n", MACHINE_LABEL);
    printf("=====================================\n");
    printf("Output path: %sTierB/\n", PATH);
    printf("Nx=%d Ny=%d dx=%.12lf dy=%.12lf dt=%.12lf iter=%d max_time=%.12lf\n",
           Nx, Ny, dx, dy, dt, iter, iter*dt);
    printf("Parameters: theta=%.12lf m=%.12lf eps=%.12lf alpha=%.12lf ell=%.12lf dA=%.12lf dF=%.12lf beta_B2=%.12lf\n",
           theta, (double)m, eps, alpha, l, dA, dF, B2_BETA);
    printf("Criteria: Mbar < %.6e, M_inf < %.6e at Tobs=%.6lf\n", MSTAR, FSTAR, B2_TOBS);

    make_output_dirs();

    auto t0 = std::chrono::steady_clock::now();
    run_B2_only();
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    std::string runtime_file = std::string(PATH) + "TierB/summary/runtime_B2_" + MACHINE_LABEL + ".txt";
    FILE* fr = fopen(runtime_file.c_str(), "w");
    if (fr) {
        fprintf(fr, "# Wall-clock time for B2-only driver\n");
        fprintf(fr, "machine %s\n", MACHINE_LABEL);
        fprintf(fr, "elapsed_seconds %.12lf\n", elapsed);
        fclose(fr);
    }

    printf("\n=====================================\n");
    printf("B2-only computations finished (%s).\n", MACHINE_LABEL);
    printf("Output directory: %sTierB/summary\n", PATH);
    printf("Elapsed wall time: %.3lf seconds\n", elapsed);
    printf("=====================================\n");
    return 0;
}
