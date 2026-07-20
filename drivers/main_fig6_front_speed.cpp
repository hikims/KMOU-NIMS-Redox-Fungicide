// =============================================================================
//  main6.cpp  --  Figure 6: recovered-state front speed
// =============================================================================
//
//  PURPOSE
//    Verify the spreading speed of Theorem 3.5,
//
//        c_inf = 2 * sqrt( dF * (1 - B) / (1 + theta^m) ),   B = beta * theta^m.
//
//    Every other analytical result in the paper already has a verification
//    panel (Proposition 3.3 -> Fig 4a, redox recovery -> Fig 4b, Theorem 3.4
//    -> Fig 4d).  This run supplies the missing one.
//
//  WHAT IS DIFFERENT FROM main_fig5.cpp
//    The physics, the discretisation, the parameters and the initial condition
//    are IDENTICAL.  Only two things change:
//
//      1. the run is longer          (t_final = 160 instead of 100)
//      2. a scalar diagnostic is recorded instead of field snapshots
//
//    The diagnostic is the number of cells whose fungal density exceeds a
//    threshold.  Because the lesion stays essentially circular about (5,5),
//    its area determines an equivalent radius,
//
//        r(t) = sqrt( n_cells(t) * dx * dy / pi ),
//
//    and the front speed is the slope of r(t).  Counting cells is free, needs
//    no interpolation, and costs one O(Nx*Ny) pass per output time.
//
//    Dumping fields at dt_out = 1.0 would produce roughly 12 GB across the
//    three doses.  This diagnostic produces about 10 kB.
//
//  WHY THREE THRESHOLDS
//    A reviewer will ask whether the measured speed depends on the level used
//    to define the front.  The three levels are 0.25, 0.50 and 0.75 times the
//    recovered plateau (1 - B), so the middle one is exactly the level
//    F_lev = (1 - B)/2 quoted in Section 5.6.  They scale automatically with
//    beta.  The actual values are written into the file header, so the
//    plotting script never has to assume them.
//
//  CELL CONVENTION
//    Identical to fungal_area_percent() in main_fig5.cpp: the domain is split
//    into (Nx-1) x (Ny-1) cells, each carrying the four-point average of its
//    corners, each of area dx*dy.  Using the same convention as the paper's
//    infected-area functional I(t) keeps Figure 6 consistent with Figure 2.
//
//  OUTPUT
//    <PATH>front_speed/front_beta_<beta>_A0_<A0>.dat
//        header  : Nx Ny dx dy levels ...   (read by scripts/plot_fig6.py)
//        columns : t  n0 n1 n2  r0 r1 r2  total_F  maxF  meanR
//
//    <PATH>front_speed/summary_fig6.txt
//        one line per dose with the final radii, for a quick sanity check.
//
//  BUILD / RUN
//    ./exec_fig6.sh                  # builds and runs all three doses
//    ./main6                         # beta = 10, A0 = 0, 2.25, 5
//    ./main6 10 0 2.25 5             # explicit
//    ./main6 8                       # secondary check at beta = 8
//
// =============================================================================

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
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

// Output root.  The Makefile passes this at compile time:
//     make main_fig6_front_speed OUTPUT_ROOT=./results/
// The value must end with a slash.  Keep this identical to the other
// drivers in drivers/ so every target writes into the same tree.
#ifndef OUTPUT_ROOT
#define OUTPUT_ROOT "./results/"
#endif

#ifndef MACHINE_LABEL
#define MACHINE_LABEL "local"
#endif

#define PATH OUTPUT_ROOT

// -----------------------------------------------------------------------------
// Figure 6 run settings
// -----------------------------------------------------------------------------
//  t_final = 160
//      A front starting at radius 2 reaches the boundary at
//      t_wall = (Lx/2 - 2) / c_inf ~ 157 for beta = 10.  Running to 160 gives
//      the full fit window [40, 140] plus a little of the confined regime.
//      Stopping at 100, as the Figure 5 run does, leaves no usable window.
//
//  dt_out = 1.0
//      101 samples inside the fit window.  Panel (b) differentiates r(t)
//      numerically, so a coarse output stride would dominate the result.
//      The diagnostic is scalar, so this costs essentially nothing.
// -----------------------------------------------------------------------------
static const double fig6_t_final = 160.0;
static const double fig6_dt_out  = 1.0;

static const double fig6_beta_default = 10.0;
static const double fig6_A0_values[]  = {0.0, 2.25, 5.0};
static const int    n_fig6_A0_values  =
    sizeof(fig6_A0_values) / sizeof(fig6_A0_values[0]);

// Front thresholds as fractions of the recovered plateau (1 - B).
static const int    n_front_levels = 3;
static const double front_level_fractions[n_front_levels] = {0.25, 0.50, 0.75};

// Initial lesion, identical to the Figure 5 run.
static const double patch_cx     = 5.0;
static const double patch_cy     = 5.0;
static const double patch_radius = 2.0;
static const double patch_inside = 0.5;


// =============================================================================
// Directory helpers (identical to main_fig5.cpp)
// =============================================================================
static void ensure_dir(const char* path) {
    struct stat st;
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0755) != 0 && errno != EEXIST) {
            perror("Fail to make directory");
            exit(1);
        }
    } else if (!S_ISDIR(st.st_mode)) {
        fprintf(stderr, "Wrong directory path: %s\n", path);
        exit(1);
    }
}

static void ensure_dir_recursive(const char* path) {
    if (path == nullptr || path[0] == '\0') return;

    char tmp[1024];
    size_t len = strlen(path);
    if (len >= sizeof(tmp)) {
        fprintf(stderr, "Directory path is too long: %s\n", path);
        exit(1);
    }

    strcpy(tmp, path);
    if (len > 1 && tmp[len - 1] == '/') tmp[len - 1] = '\0';

    for (char* p = tmp + 1; *p; ++p) {
        if (*p == '/') {
            *p = '\0';
            ensure_dir(tmp);
            *p = '/';
        }
    }
    ensure_dir(tmp);
}


// =============================================================================
// Observables
// =============================================================================
static double compute_total_mass(double **U) {
    int a = 0, b = Nx - 1, c = 0, d = Ny - 1;
    double first_term = 0.0, second_term = 0.0;
    double third_term = 0.0, fourth_term = 0.0, fifth_term = 0.0;

    for (int j = 1; j < Ny - 1; j++) {
        first_term  += U[a][j];
        second_term += U[b][j];
    }
    for (int i = 1; i < Nx - 1; i++) {
        third_term  += U[i][c];
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

// Number of cells whose four-point average exceeds each threshold.
// Same cell convention as fungal_area_percent() in main_fig5.cpp.
static void count_super_level_cells(double **F,
                                    const double *levels,
                                    long *counts) {
    for (int q = 0; q < n_front_levels; q++) counts[q] = 0;

    for (int i = 0; i < Nx - 1; i++) {
        for (int j = 0; j < Ny - 1; j++) {
            const double Fav = 0.25 * (F[i][j] + F[i+1][j]
                                       + F[i][j+1] + F[i+1][j+1]);
            for (int q = 0; q < n_front_levels; q++) {
                if (Fav >= levels[q]) counts[q]++;
            }
        }
    }
}

static double equivalent_radius(long n_cells) {
    return sqrt((double)n_cells * dx * dy / M_PI);
}

static double field_max(double **U) {
    double v = U[0][0];
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            if (U[i][j] > v) v = U[i][j];
    return v;
}

static double field_mean(double **U) {
    double s = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            s += U[i][j];
    return s / ((double)Nx * (double)Ny);
}


// =============================================================================
// Initial condition (identical to main_fig5.cpp)
// =============================================================================
static void initialize_fungal_patch(double **F) {
    const double r2 = patch_radius * patch_radius;

    for (int i = 0; i < Nx; i++) {
        const double x = i * dx;
        for (int j = 0; j < Ny; j++) {
            const double y = j * dy;
            const double dist2 = (x - patch_cx) * (x - patch_cx)
                               + (y - patch_cy) * (y - patch_cy);
            F[i][j] = (dist2 <= r2) ? patch_inside : 0.0;
        }
    }
}


// =============================================================================
// One simulation
// =============================================================================
static void run_simulation(double A0, double beta_val, FILE *fSummary) {
    // ---------- derived analytic quantities ----------
    const double B      = beta_val * pow(theta, m);
    const double a_inf  = (1.0 - B) / (1.0 + pow(theta, m));
    const double c_inf  = (a_inf > 0.0) ? 2.0 * sqrt(dF * a_inf) : 0.0;
    const double width  = (a_inf > 0.0) ? sqrt(dF / a_inf) : 0.0;
    const double Lx     = (Nx - 1) * dx;
    const double Ly     = (Ny - 1) * dy;
    const double t_wall = (c_inf > 0.0)
                        ? (0.5 * Lx - patch_radius) / c_inf : 0.0;

    double levels[n_front_levels];
    for (int q = 0; q < n_front_levels; q++) {
        levels[q] = front_level_fractions[q] * (1.0 - B);
    }

    // ---------- memory ----------
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

    // ---------- ADI diffusion coefficients ----------
    double rAx = (dA * 0.5 * dt) / (dx * dx);
    double rAy = (dA * 0.5 * dt) / (dy * dy);
    double rFx = (dF * 0.5 * dt) / (dx * dx);
    double rFy = (dF * 0.5 * dt) / (dy * dy);

    // ---------- initial condition ----------
    srand(1);
    initialization(R, Nx, Ny, initR);
    initialization(A, Nx, Ny, A0);
    initialize_fungal_patch(F);

    // ---------- LU setup (Neumann BC) ----------
    for (int i = 0; i < Nx; i++) {
        A_Ax[i] = 1 + 2 * rAx;  B_Ax[i] = -rAx;  C_Ax[i] = -rAx;
    }
    A_Ax[0] = 1 + 2 * rAx;       B_Ax[0] = 0.0;         C_Ax[0] = -2 * rAx;
    A_Ax[Nx - 1] = 1 + 2 * rAx;  B_Ax[Nx - 1] = -2 * rAx; C_Ax[Nx - 1] = 0.0;

    for (int j = 0; j < Ny; j++) {
        A_Ay[j] = 1 + 2 * rAy;  B_Ay[j] = -rAy;  C_Ay[j] = -rAy;
    }
    A_Ay[0] = 1 + 2 * rAy;       B_Ay[0] = 0.0;         C_Ay[0] = -2 * rAy;
    A_Ay[Ny - 1] = 1 + 2 * rAy;  B_Ay[Ny - 1] = -2 * rAy; C_Ay[Ny - 1] = 0.0;

    LUx(A_Ax, B_Ax, C_Ax, L_Ax, U_Ax, Nx);
    LUy(A_Ay, B_Ay, C_Ay, L_Ay, U_Ay, Ny);

    for (int i = 0; i < Nx; i++) {
        A_Fx[i] = 1 + 2 * rFx;  B_Fx[i] = -rFx;  C_Fx[i] = -rFx;
    }
    A_Fx[0] = 1 + 2 * rFx;       B_Fx[0] = 0.0;         C_Fx[0] = -2 * rFx;
    A_Fx[Nx - 1] = 1 + 2 * rFx;  B_Fx[Nx - 1] = -2 * rFx; C_Fx[Nx - 1] = 0.0;

    for (int j = 0; j < Ny; j++) {
        A_Fy[j] = 1 + 2 * rFy;  B_Fy[j] = -rFy;  C_Fy[j] = -rFy;
    }
    A_Fy[0] = 1 + 2 * rFy;       B_Fy[0] = 0.0;         C_Fy[0] = -2 * rFy;
    A_Fy[Ny - 1] = 1 + 2 * rFy;  B_Fy[Ny - 1] = -2 * rFy; C_Fy[Ny - 1] = 0.0;

    LUx(A_Fx, B_Fx, C_Fx, L_Fx, U_Fx, Nx);
    LUy(A_Fy, B_Fy, C_Fy, L_Fy, U_Fy, Ny);

    // ---------- output file ----------
    ensure_dir_recursive(PATH);
    ensure_dir_recursive(PATH "front_speed");

    char front_file[1024];
    snprintf(front_file, sizeof(front_file),
             PATH "front_speed/front_beta_%.3f_A0_%.3f.dat", beta_val, A0);

    FILE *fFront = fopen(front_file, "w");
    if (!fFront) {
        perror("fail to open front-speed file");
        exit(1);
    }

    // The header carries the grid metadata and the thresholds, so the
    // plotting script never has to hard-code dx or the levels.
    fprintf(fFront, "# Figure 6 front-speed diagnostic\n");
    fprintf(fFront, "# A0     = %.12f\n", A0);
    fprintf(fFront, "# beta   = %.12f\n", beta_val);
    fprintf(fFront, "# B      = %.12f\n", B);
    fprintf(fFront, "# a_inf  = %.12f\n", a_inf);
    fprintf(fFront, "# c_inf  = %.12f\n", c_inf);
    fprintf(fFront, "# width  = %.12f   # sqrt(dF/a_inf)\n", width);
    fprintf(fFront, "# t_wall = %.12f\n", t_wall);
    fprintf(fFront, "# Nx = %d\n# Ny = %d\n", Nx, Ny);
    fprintf(fFront, "# dx = %.12f\n# dy = %.12f\n", dx, dy);
    fprintf(fFront, "# Lx = %.12f\n# Ly = %.12f\n", Lx, Ly);
    fprintf(fFront, "# dt = %.12f\n# dt_out = %.12f\n", dt, fig6_dt_out);
    fprintf(fFront, "# t_final = %.12f\n", fig6_t_final);
    fprintf(fFront, "# levels =");
    for (int q = 0; q < n_front_levels; q++) fprintf(fFront, " %.12f", levels[q]);
    fprintf(fFront, "\n");
    fprintf(fFront, "# primary_level_index = 1   # 0.5*(1-B), Section 5.6\n");
    fprintf(fFront, "# cell_convention = four-point average over (Nx-1)x(Ny-1) cells\n");
    fprintf(fFront, "# columns: t n0 n1 n2 r0 r1 r2 total_F maxF meanR\n");
    fflush(fFront);

    // ---------- stepping ----------
    const long n_steps  = llround(fig6_t_final / dt);
    const long out_every = llround(fig6_dt_out / dt);

    if (out_every < 1) {
        fprintf(stderr, "dt_out = %.6f is smaller than dt = %.6f\n",
                fig6_dt_out, dt);
        exit(1);
    }

    long counts[n_front_levels];
    bool wall_warned = false;

    // Record t = 0 before any step is taken.
    {
        count_super_level_cells(F, levels, counts);
        fprintf(fFront, "%.6f", 0.0);
        for (int q = 0; q < n_front_levels; q++) fprintf(fFront, " %ld", counts[q]);
        for (int q = 0; q < n_front_levels; q++)
            fprintf(fFront, " %.12f", equivalent_radius(counts[q]));
        fprintf(fFront, " %.12f %.12f %.12f\n",
                compute_total_mass(F), field_max(F), field_mean(R));
        fflush(fFront);
    }

    tbb::task_arena arena;

    for (long k = 0; k < n_steps; k++)
    {
        // Neumann BC
        for (int j = 0; j < Ny; j++) {
            A[0][j]      = A[1][j];
            A[Nx - 1][j] = A[Nx - 2][j];
            F[0][j]      = F[1][j];
            F[Nx - 1][j] = F[Nx - 2][j];
        }
        for (int i = 0; i < Nx; i++) {
            A[i][0]      = A[i][1];
            A[i][Ny - 1] = A[i][Ny - 2];
            F[i][0]      = F[i][1];
            F[i][Ny - 1] = F[i][Ny - 2];
        }

        // x-direction
        arena.execute([&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, Ny),
                              [&](const tbb::blocked_range<int> &r) {
                for (int j = r.begin(); j < r.end(); j++) {
                    for (int i = 0; i < Nx; i++) {
                        int jp = (j == Ny - 1) ? j - 1 : j + 1;
                        int jm = (j == 0)      ? j + 1 : j - 1;

                        Ax[j][i] = coefficient_dA(A[i][j], rAy,
                                                  A[i][jp], A[i][j], A[i][jm],
                                                  F[i][j], alpha);

                        // Runtime-beta overload, as in main4.cpp
                        Fx[j][i] = coefficient_dF(F[i][j], rFy,
                                                  F[i][jp], F[i][j], F[i][jm],
                                                  R[i][j], m, theta, beta_val);
                    }

                    solx(L_Ax, U_Ax, &Ax[j][0], &z_Ax[j][0], &x_Ax[j][0], Nx);
                    solx(L_Fx, U_Fx, &Fx[j][0], &z_Fx[j][0], &x_Fx[j][0], Nx);

                    for (int i = 0; i < Nx; i++) {
                        A_new[i][j] = x_Ax[j][i];
                        F_new[i][j] = x_Fx[j][i];
                        R_new[i][j] = R[i][j]
                                    + dt * (eps * (1.0 - R[i][j]) - A[i][j] * R[i][j]);
                    }
                }
            });
        });

        // y-direction
        arena.execute([&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, Nx),
                              [&](const tbb::blocked_range<int> &r) {
                for (int i = r.begin(); i < r.end(); i++) {
                    for (int j = 0; j < Ny; j++) {
                        int ip = (i == Nx - 1) ? i - 1 : i + 1;
                        int im = (i == 0)      ? i + 1 : i - 1;

                        Ay[i][j] = coefficient_dA(A_new[i][j], rAx,
                                                  A_new[ip][j], A_new[i][j], A_new[im][j],
                                                  F_new[i][j], alpha);

                        Fy[i][j] = coefficient_dF(F_new[i][j], rFx,
                                                  F_new[ip][j], F_new[i][j], F_new[im][j],
                                                  R_new[i][j], m, theta, beta_val);
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

        const long step = k + 1;

        if (step % out_every == 0) {
            const double current_time = (double)step * dt;

            count_super_level_cells(F, levels, counts);

            const double total_F = compute_total_mass(F);
            const double maxF    = field_max(F);
            const double meanR   = field_mean(R);

            fprintf(fFront, "%.6f", current_time);
            for (int q = 0; q < n_front_levels; q++) fprintf(fFront, " %ld", counts[q]);
            for (int q = 0; q < n_front_levels; q++)
                fprintf(fFront, " %.12f", equivalent_radius(counts[q]));
            fprintf(fFront, " %.12f %.12f %.12f\n", total_F, maxF, meanR);
            fflush(fFront);

            const double r_mid = equivalent_radius(counts[1]);

            printf("[A0=%.3f beta=%.3f] t=%7.2f  r=%.5f  n=%ld  maxF=%.5f\n",
                   A0, beta_val, current_time, r_mid, counts[1], maxF);

            // The KPP speed is a free-space statement.  Once the lesion feels
            // the boundary, the slope is no longer c_inf and must not be used
            // in the fit.
            if (!wall_warned && r_mid > 0.5 * Lx - 5.0 * dx) {
                printf("  NOTE: the front is within five cells of the boundary "
                       "at t = %.2f. Exclude later times from the fit.\n",
                       current_time);
                wall_warned = true;
            }

            if (isnan(maxF) || maxF > 1.0 + 1e-9) {
                printf("  Unstable at t = %.6f: maxF = %.6f\n", current_time, maxF);
                break;
            }
        }
    }

    fclose(fFront);

    // ---------- per-dose summary line ----------
    count_super_level_cells(F, levels, counts);
    fprintf(fSummary, "%.6f %.6f %.12f %.12f %.12f %.12f\n",
            A0, beta_val, c_inf, t_wall,
            equivalent_radius(counts[1]), compute_total_mass(F));
    fflush(fSummary);

    printf("\nFinished A0 = %.3f\n", A0);
    printf("  final r (level %.4f) = %.5f\n", levels[1], equivalent_radius(counts[1]));
    printf("  written to %s\n", front_file);

    // ---------- free ----------
    freeMallocDouble(L_Ax, Lp);
    freeMallocDouble(L_Ay, Lp);
    freeMallocDouble(L_Fx, Lp);
    freeMallocDouble(L_Fy, Lp);

    freeMallocDouble(A, Nx);
    freeMallocDouble(A_new, Nx);
    freeMallocDouble(R, Nx);
    freeMallocDouble(R_new, Nx);
    freeMallocDouble(F, Nx);
    freeMallocDouble(F_new, Nx);

    freeMallocDouble(x_Ax, Ny);
    freeMallocDouble(x_Ay, Nx);
    freeMallocDouble(z_Ax, Ny);
    freeMallocDouble(z_Ay, Nx);
    freeMallocSingle(U_Ax);
    freeMallocSingle(U_Ay);
    freeMallocSingle(A_Ax);
    freeMallocSingle(B_Ax);
    freeMallocSingle(C_Ax);
    freeMallocSingle(A_Ay);
    freeMallocSingle(B_Ay);
    freeMallocSingle(C_Ay);

    freeMallocDouble(x_Fx, Ny);
    freeMallocDouble(x_Fy, Nx);
    freeMallocDouble(z_Fx, Ny);
    freeMallocDouble(z_Fy, Nx);
    freeMallocSingle(U_Fx);
    freeMallocSingle(U_Fy);
    freeMallocSingle(A_Fx);
    freeMallocSingle(B_Fx);
    freeMallocSingle(C_Fx);
    freeMallocSingle(A_Fy);
    freeMallocSingle(B_Fy);
    freeMallocSingle(C_Fy);

    freeMallocDouble(Ax, Ny);
    freeMallocDouble(Ay, Nx);
    freeMallocDouble(Fx, Ny);
    freeMallocDouble(Fy, Nx);
}


// =============================================================================
// main
// =============================================================================
int main(int argc, char **argv)
{
    double beta_val = fig6_beta_default;
    std::vector<double> doses;

    if (argc >= 2) beta_val = atof(argv[1]);
    for (int i = 2; i < argc; i++) doses.push_back(atof(argv[i]));
    if (doses.empty()) {
        for (int i = 0; i < n_fig6_A0_values; i++) doses.push_back(fig6_A0_values[i]);
    }

    const double B      = beta_val * pow(theta, m);
    const double a_inf  = (1.0 - B) / (1.0 + pow(theta, m));
    const double Lx     = (Nx - 1) * dx;

    printf("\n=====================================\n");
    printf("Figure 6: recovered-state front speed\n");
    printf("=====================================\n");
    printf("  beta        = %.6f\n", beta_val);
    printf("  B           = %.6f\n", B);

    if (B >= 1.0) {
        fprintf(stderr,
                "\nERROR: B = %.6f >= 1. The recovered state is extinction and\n"
                "no propagating front exists. Use beta < %.6f.\n",
                B, 1.0 / pow(theta, m));
        return 1;
    }

    const double c_inf  = 2.0 * sqrt(dF * a_inf);
    const double width  = sqrt(dF / a_inf);
    const double t_wall = (0.5 * Lx - patch_radius) / c_inf;

    printf("  a_inf       = %.6f\n", a_inf);
    printf("  c_inf       = %.6f\n", c_inf);
    printf("  front width = %.6f  (%.1f cells)\n", width, width / dx);
    printf("  t_wall      = %.1f\n", t_wall);
    printf("  t_final     = %.1f,  dt_out = %.2f\n", fig6_t_final, fig6_dt_out);
    printf("  grid        = %d x %d,  dx = %.6f,  Lx = %.4f\n", Nx, Ny, dx, Lx);
    printf("  output root = %s\n", PATH);

    if (width / dx < 4.0) {
        printf("\n  WARNING: the front spans only %.1f cells. The speed is not\n"
               "  resolved at this dx; expect a biased slope.\n", width / dx);
    }
    if (t_wall < 0.9 * fig6_t_final) {
        printf("\n  WARNING: the front reaches the boundary at t = %.1f, before\n"
               "  t_final = %.1f. Restrict the fit to earlier times.\n",
               t_wall, fig6_t_final);
    }

    printf("=====================================\n");

    ensure_dir_recursive(PATH);
    ensure_dir_recursive(PATH "front_speed");

    char summary_file[1024];
    snprintf(summary_file, sizeof(summary_file),
             PATH "front_speed/summary_fig6.txt");

    FILE *fSummary = fopen(summary_file, "w");
    if (!fSummary) {
        perror("fail to open summary file");
        return 1;
    }
    fprintf(fSummary, "# A0 beta c_inf t_wall final_r total_F\n");
    fflush(fSummary);

    for (double A0 : doses) {
        printf("\n-------------------------------------\n");
        printf("Running: A0 = %.6f, beta = %.6f\n", A0, beta_val);
        printf("-------------------------------------\n");
        run_simulation(A0, beta_val, fSummary);
    }

    fclose(fSummary);

    printf("\nAll Figure 6 runs finished.\n");
    printf("Now plot with:\n");
    printf("  python scripts/plot_fig6.py --beta %g\n", beta_val);
    printf("Then put the fitted speed into the two \\todo{} marks in\n");
    printf("Section 5.6 of main.tex and delete the macro.\n");

    return 0;
}
