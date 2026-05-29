/*
 * Driver for main-manuscript Figures 3 and 4.
 * Generates time-series data for phase-transition dynamics and asymptotic-threshold verification.
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
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>

#ifndef OUTPUT_ROOT
#define OUTPUT_ROOT "./results/"
#endif
#define PATH             OUTPUT_ROOT
#define TIME_SERIES_PATH PATH "time_series/"

// =============================================================
//  Usage:
//    ./main4 <A0>            -> uses default beta from parameters.h
//    ./main4 <A0> <beta>     -> overrides beta at runtime
// =============================================================


// -------------------------------------------------------------
// Directory helpers
// -------------------------------------------------------------
static void ensure_dir(const char* path) {
    struct stat st;
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0755) != 0 && errno != EEXIST) {
            perror("Fail to make directory");
            fprintf(stderr, "  path = %s\n", path);
            exit(1);
        }
    } else if (!S_ISDIR(st.st_mode)) {
        fprintf(stderr, "Wrong directory path: %s\n", path);
        exit(1);
    }
}

static void ensure_dir_recursive(const char* path) {
    char tmp[1024];
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (len > 0 && tmp[len - 1] == '/') tmp[len - 1] = '\0';
    for (char* p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            ensure_dir(tmp);
            *p = '/';
        }
    }
    ensure_dir(tmp);
}


// -------------------------------------------------------------
// Field statistics helpers
// -------------------------------------------------------------
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
    double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0;
    for (int j = 1; j < Ny - 1; j++) { s1 += U[a][j]; s2 += U[b][j]; }
    for (int i = 1; i < Nx - 1; i++) { s3 += U[i][c]; s4 += U[i][d]; }
    for (int i = 1; i < Nx - 1; i++)
        for (int j = 1; j < Ny - 1; j++) s5 += U[i][j];
    return 0.25 * dx * dy *
           (U[a][c] + U[a][d] + U[b][c] + U[b][d]
            + 2.0 * (s1 + s2 + s3 + s4)
            + 4.0 * s5);
}

static double compute_mean_value(double **U) {
    double sum = 0.0;
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++) sum += U[i][j];
    return sum / (double)(Nx * Ny);
}

static double compute_min_value(double **U) {
    double vmin = U[0][0];
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            if (U[i][j] < vmin) vmin = U[i][j];
    return vmin;
}

static double compute_max_value(double **U) {
    double vmax = U[0][0];
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            if (U[i][j] > vmax) vmax = U[i][j];
    return vmax;
}


// -------------------------------------------------------------
// Snapshot writer
// -------------------------------------------------------------
static void save_field_snapshot(const char* filepath,
                                double **A, double **F, double **R) {
    FILE *fp = fopen(filepath, "w");
    if (!fp) { perror("fail to open snapshot file"); exit(1); }
    fprintf(fp, "# x y A F R\n");
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            fprintf(fp, "%lf %lf %.12lf %.12lf %.12lf\n",
                    i * dx, j * dy, A[i][j], F[i][j], R[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


// -------------------------------------------------------------
// Theorem 3.8 prediction
// -------------------------------------------------------------
static void compute_theory(double beta_val,
                           double *F_asymptote, int *case_id,
                           double *f1_out, double *g1_out)
{
    double theta_m = std::pow(theta, m);
    double f1 = 1.0 / (1.0 + theta_m);
    double g1 = beta_val * theta_m / (1.0 + theta_m);
    if (f1_out) *f1_out = f1;
    if (g1_out) *g1_out = g1;

    if (f1 - g1 < 0.0) {
        *case_id = 1;            // F -> 0
        *F_asymptote = 0.0;
    } else {
        *case_id = 2;            // F -> 1 - g(1)/f(1)
        *F_asymptote = 1.0 - g1 / f1;
    }
}


// -------------------------------------------------------------
// Main simulation
// -------------------------------------------------------------
static double run_simulation(double A0, double beta_val, FILE* fGlobal)
{
    // ---------- Memory allocation ----------
    double **A     = mallocDouble(Nx, Ny);
    double **R     = mallocDouble(Nx, Ny);
    double **F     = mallocDouble(Nx, Ny);
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
    double *U_Ax  = mallocSingle(Nx);
    double *U_Ay  = mallocSingle(Ny);

    double **L_Fx = mallocDouble(Lp, Nx);
    double **L_Fy = mallocDouble(Lp, Ny);
    double *U_Fx  = mallocSingle(Nx);
    double *U_Fy  = mallocSingle(Ny);

    // ---------- ADI diffusion coefficients ----------
    double rAx = (dA * 0.5 * dt) / (dx * dx);
    double rAy = (dA * 0.5 * dt) / (dy * dy);
    double rFx = (dF * 0.5 * dt) / (dx * dx);
    double rFy = (dF * 0.5 * dt) / (dy * dy);

    // ---------- Initial conditions ----------
    srand(1);
    initialization(R, Nx, Ny, initR);
    initialization(A, Nx, Ny, A0);
    Circular(F, Nx, Ny, dx, dy,
             F_center_cx, F_center_cy, F_radius,
             F_inside, F_outside);

    // ---------- LU decomposition ----------
    for (int i = 0; i < Nx; i++) { A_Ax[i] = 1 + 2*rAx; B_Ax[i] = -rAx; C_Ax[i] = -rAx; }
    A_Ax[0] = 1 + 2*rAx; B_Ax[0] = 0.0;     C_Ax[0] = -2*rAx;
    A_Ax[Nx-1] = 1 + 2*rAx; B_Ax[Nx-1] = -2*rAx; C_Ax[Nx-1] = 0.0;
    for (int j = 0; j < Ny; j++) { A_Ay[j] = 1 + 2*rAy; B_Ay[j] = -rAy; C_Ay[j] = -rAy; }
    A_Ay[0] = 1 + 2*rAy; B_Ay[0] = 0.0;     C_Ay[0] = -2*rAy;
    A_Ay[Ny-1] = 1 + 2*rAy; B_Ay[Ny-1] = -2*rAy; C_Ay[Ny-1] = 0.0;
    LUx(A_Ax, B_Ax, C_Ax, L_Ax, U_Ax, Nx);
    LUy(A_Ay, B_Ay, C_Ay, L_Ay, U_Ay, Ny);

    for (int i = 0; i < Nx; i++) { A_Fx[i] = 1 + 2*rFx; B_Fx[i] = -rFx; C_Fx[i] = -rFx; }
    A_Fx[0] = 1 + 2*rFx; B_Fx[0] = 0.0;     C_Fx[0] = -2*rFx;
    A_Fx[Nx-1] = 1 + 2*rFx; B_Fx[Nx-1] = -2*rFx; C_Fx[Nx-1] = 0.0;
    for (int j = 0; j < Ny; j++) { A_Fy[j] = 1 + 2*rFy; B_Fy[j] = -rFy; C_Fy[j] = -rFy; }
    A_Fy[0] = 1 + 2*rFy; B_Fy[0] = 0.0;     C_Fy[0] = -2*rFy;
    A_Fy[Ny-1] = 1 + 2*rFy; B_Fy[Ny-1] = -2*rFy; C_Fy[Ny-1] = 0.0;
    LUx(A_Fx, B_Fx, C_Fx, L_Fx, U_Fx, Nx);
    LUy(A_Fy, B_Fy, C_Fy, L_Fy, U_Fy, Ny);

    // ---------- Output folders ----------
    ensure_dir_recursive(PATH);
    ensure_dir_recursive(TIME_SERIES_PATH);
    ensure_dir_recursive(TIME_SERIES_PATH "A_R_F_timeseries");
    ensure_dir_recursive(TIME_SERIES_PATH "spatial_fields");
    ensure_dir_recursive(TIME_SERIES_PATH "theorem_3_8");

    char run_tag[128];
    snprintf(run_tag, sizeof(run_tag), "A0_%.3f_beta_%.3f", A0, beta_val);

    char dir_spatial[1024];
    snprintf(dir_spatial, sizeof(dir_spatial),
             TIME_SERIES_PATH "spatial_fields/%s", run_tag);
    ensure_dir(dir_spatial);

    char ts_file[1024];
    snprintf(ts_file, sizeof(ts_file),
             TIME_SERIES_PATH "A_R_F_timeseries/%s.dat", run_tag);
    FILE *fts = fopen(ts_file, "w");
    if (!fts) { perror("fail to open ts_file"); exit(1); }
    fprintf(fts, "# Time series of A, R, F (mean / min / max over Omega)\n");
    fprintf(fts, "# A0 = %.6f, beta = %.6f\n", A0, beta_val);
    fprintf(fts, "# col1=t  col2=mean_A  col3=mean_R  col4=mean_F  "
                 "col5=min_A col6=max_A  col7=min_R col8=max_R  "
                 "col9=min_F col10=max_F  col11=fungal_area_percent  "
                 "col12=total_A col13=total_F\n");

    double F_asymptote, f1, g1;
    int case_id;
    compute_theory(beta_val, &F_asymptote, &case_id, &f1, &g1);

    char th_evol_file[1024];
    snprintf(th_evol_file, sizeof(th_evol_file),
             TIME_SERIES_PATH "theorem_3_8/%s_F_evolution.dat", run_tag);
    FILE *fth = fopen(th_evol_file, "w");
    if (!fth) { perror("fail to open th_evol_file"); exit(1); }
    fprintf(fth, "# Theorem 3.8 validation : F evolution vs asymptote\n");
    fprintf(fth, "# A0 = %.6f, beta = %.6f, theta = %.6f, m = %.6f\n",
            A0, beta_val, theta, m);
    fprintf(fth, "# beta*theta^m = %.12f\n", beta_val * std::pow(theta, m));
    fprintf(fth, "# Case = %d  (1: F->0,  2: F->1-g(1)/f(1))\n", case_id);
    fprintf(fth, "# F_asymptote = %.12f\n", F_asymptote);
    fprintf(fth, "# col1=t  col2=mean_F  col3=max_F  col4=F_asymptote  col5=|mean_F-F_asymptote|\n");

    int eval_steps[n_eval_times];
    bool eval_done[n_eval_times];
    for (int e = 0; e < n_eval_times; e++) {
        eval_steps[e] = static_cast<int>(eval_times[e] / dt);
        eval_done[e]  = false;
    }

    tbb::task_arena arena;

    // ---------- Time iteration ----------
    for (int k = 0; k < iter; k++)
    {
        // Neumann BC
        for (int j = 0; j < Ny; j++) {
            A[0][j] = A[1][j];     A[Nx-1][j] = A[Nx-2][j];
            F[0][j] = F[1][j];     F[Nx-1][j] = F[Nx-2][j];
        }
        for (int i = 0; i < Nx; i++) {
            A[i][0] = A[i][1];     A[i][Ny-1] = A[i][Ny-2];
            F[i][0] = F[i][1];     F[i][Ny-1] = F[i][Ny-2];
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
                        // <-- using overloaded coefficient_dF with runtime beta
                        Fx[j][i] = coefficient_dF(F[i][j], rFy,
                                                  F[i][jp], F[i][j], F[i][jm],
                                                  R[i][j], m, theta, beta_val);
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
                        // <-- using overloaded coefficient_dF with runtime beta
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

        if (k % cut == 0) {
            double t       = k * dt;
            double mean_A  = compute_mean_value(A);
            double mean_R  = compute_mean_value(R);
            double mean_F  = compute_mean_value(F);
            double min_A   = compute_min_value(A);
            double max_A   = compute_max_value(A);
            double min_R   = compute_min_value(R);
            double max_R   = compute_max_value(R);
            double min_F   = compute_min_value(F);
            double max_F   = compute_max_value(F);
            double area    = fungal_area_percent(F);
            double total_A = compute_total_mass(A);
            double total_F = compute_total_mass(F);

            fprintf(fts,
                "%.6lf %.12lf %.12lf %.12lf "
                "%.12lf %.12lf %.12lf %.12lf "
                "%.12lf %.12lf %.6lf "
                "%.12lf %.12lf\n",
                t, mean_A, mean_R, mean_F,
                min_A, max_A, min_R, max_R,
                min_F, max_F, area,
                total_A, total_F);
            fflush(fts);

            fprintf(fth, "%.6lf %.12lf %.12lf %.12lf %.12lf\n",
                    t, mean_F, max_F, F_asymptote, fabs(mean_F - F_asymptote));
            fflush(fth);

            printf("[A0=%.3f beta=%.3f] t=%8.3f  meanA=%.5f  meanR=%.5f  meanF=%.5f  area%%=%6.2f\n",
                   A0, beta_val, t, mean_A, mean_R, mean_F, area);

            if (F[Nx / 2][Ny / 2] > 1.0 || isnan(F[Nx / 2][Ny / 2])) {
                printf("Unstable at step %d (t=%.6lf): F=%.6lf\n",
                       k, t, F[Nx / 2][Ny / 2]);
                break;
            }
        }

        for (int e = 0; e < n_eval_times; e++) {
            if (!eval_done[e] && k >= eval_steps[e]) {
                char snap_file[1024];
                snprintf(snap_file, sizeof(snap_file),
                         "%s/field_t%05.1f.dat", dir_spatial, eval_times[e]);
                save_field_snapshot(snap_file, A, F, R);
                eval_done[e] = true;
                printf("  -> saved spatial snapshot: %s\n", snap_file);
            }
        }
    }

    // ---------- Final field & summary ----------
    char final_field_file[1024];
    snprintf(final_field_file, sizeof(final_field_file),
             "%s/field_final.dat", dir_spatial);
    save_field_snapshot(final_field_file, A, F, R);

    double final_mean_F = compute_mean_value(F);
    double final_mean_A = compute_mean_value(A);
    double final_mean_R = compute_mean_value(R);
    double final_max_F  = compute_max_value(F);
    double final_min_F  = compute_min_value(F);

    char summary_file[1024];
    snprintf(summary_file, sizeof(summary_file),
             TIME_SERIES_PATH "theorem_3_8/%s_summary.txt", run_tag);
    FILE *fsum = fopen(summary_file, "w");
    if (fsum) {
        fprintf(fsum, "Theorem 3.8 - Numerical validation\n");
        fprintf(fsum, "=====================================\n\n");
        fprintf(fsum, "Parameters\n");
        fprintf(fsum, "  A0    = %.6f\n", A0);
        fprintf(fsum, "  dA    = %.6f\n", dA);
        fprintf(fsum, "  dF    = %.6f\n", dF);
        fprintf(fsum, "  alpha = %.6f\n", alpha);
        fprintf(fsum, "  ell   = %.6f\n", l);
        fprintf(fsum, "  eps   = %.6f\n", eps);
        fprintf(fsum, "  beta  = %.6f  (runtime)\n", beta_val);
        fprintf(fsum, "  theta = %.6f\n", theta);
        fprintf(fsum, "  m     = %.6f\n\n", m);

        fprintf(fsum, "Theory\n");
        fprintf(fsum, "  beta*theta^m = %.12f\n", beta_val * std::pow(theta, m));
        fprintf(fsum, "  f(1) = %.12f\n", f1);
        fprintf(fsum, "  g(1) = %.12f\n", g1);
        fprintf(fsum, "  f(1) - g(1) = %.12f\n", f1 - g1);
        if (case_id == 1) fprintf(fsum, "  -> Case I  : F -> 0\n");
        else              fprintf(fsum, "  -> Case II : F -> 1 - g(1)/f(1)\n");
        fprintf(fsum, "  Theoretical F_asymptote = %.12f\n\n", F_asymptote);

        fprintf(fsum, "Final numerical values (t = %.4f)\n", (iter - 1) * dt);
        fprintf(fsum, "  mean_A = %.12f\n", final_mean_A);
        fprintf(fsum, "  mean_R = %.12f\n", final_mean_R);
        fprintf(fsum, "  mean_F = %.12f\n", final_mean_F);
        fprintf(fsum, "  min_F  = %.12f\n", final_min_F);
        fprintf(fsum, "  max_F  = %.12f\n", final_max_F);
        fprintf(fsum, "  |mean_F - F_asymptote| = %.12f\n",
                fabs(final_mean_F - F_asymptote));
        fclose(fsum);
    }

    if (fGlobal) {
        fprintf(fGlobal,
                "%.6f %.6f %d %.12f %.12f %.12f %.12f %.12f\n",
                A0, beta_val, case_id,
                beta_val * std::pow(theta, m),
                F_asymptote,
                final_mean_F,
                final_max_F,
                fabs(final_mean_F - F_asymptote));
        fflush(fGlobal);
    }

    fclose(fts);
    fclose(fth);

    // ---------- Free memory ----------
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
    freeMallocSingle(A_Ax); freeMallocSingle(B_Ax); freeMallocSingle(C_Ax);
    freeMallocSingle(A_Ay); freeMallocSingle(B_Ay); freeMallocSingle(C_Ay);

    freeMallocDouble(x_Fx, Ny);
    freeMallocDouble(x_Fy, Nx);
    freeMallocDouble(z_Fx, Ny);
    freeMallocDouble(z_Fy, Nx);
    freeMallocSingle(U_Fx);
    freeMallocSingle(U_Fy);
    freeMallocSingle(A_Fx); freeMallocSingle(B_Fx); freeMallocSingle(C_Fx);
    freeMallocSingle(A_Fy); freeMallocSingle(B_Fy); freeMallocSingle(C_Fy);

    freeMallocDouble(Ax, Ny);
    freeMallocDouble(Ay, Nx);
    freeMallocDouble(Fx, Ny);
    freeMallocDouble(Fy, Nx);

    return final_mean_F;
}


// -------------------------------------------------------------
// main
// -------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s A0 [beta]\n", argv[0]);
        fprintf(stderr, "  A0   : initial fungicide concentration (required)\n");
        fprintf(stderr, "  beta : optional, overrides beta from parameters.h (default %.3f)\n", beta);
        return 1;
    }

    double A0       = atof(argv[1]);
    double beta_val = (argc >= 3) ? atof(argv[2]) : beta;

    printf("\n=====================================\n");
    printf("Running simulation\n");
    printf("  A0   = %.6f\n", A0);
    printf("  beta = %.6f  %s\n", beta_val,
           (argc >= 3) ? "(from CLI)" : "(default from parameters.h)");
    printf("  beta * theta^m = %.6f\n", beta_val * std::pow(theta, m));
    printf("Output root  : %s\n", PATH);
    printf("=====================================\n");

    ensure_dir_recursive(PATH);
    ensure_dir_recursive(TIME_SERIES_PATH);
    ensure_dir_recursive(TIME_SERIES_PATH "theorem_3_8");

    char global_file[1024];
    snprintf(global_file, sizeof(global_file),
             TIME_SERIES_PATH "theorem_3_8/global_summary.txt");

    FILE *fGlobal = fopen(global_file, "a");
    if (!fGlobal) { perror("fail to open global_summary.txt"); return 1; }

    fseek(fGlobal, 0, SEEK_END);
    if (ftell(fGlobal) == 0) {
        fprintf(fGlobal,
                "# col1=A0  col2=beta  col3=case_id  col4=beta*theta^m  "
                "col5=F_asymptote  col6=final_mean_F  col7=final_max_F  "
                "col8=|final_mean_F - F_asymptote|\n");
        fflush(fGlobal);
    }

    run_simulation(A0, beta_val, fGlobal);

    fclose(fGlobal);

    printf("\nSimulation finished for A0 = %.6f, beta = %.6f\n", A0, beta_val);
    return 0;
}
