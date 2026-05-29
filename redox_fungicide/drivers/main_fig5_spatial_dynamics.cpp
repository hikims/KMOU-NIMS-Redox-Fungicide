/*
 * Driver for main-manuscript Figure 5.
 * Generates spatial snapshots at t = 0, 10, 50, and 100 for A0 = 0, 2.25, and 5.
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

// Output root for the Figure 5 rerun.
// By default, files are written under ./results/fig5_rerun_A0_0_2p25_5/.
// To change the parent directory, compile with e.g.
//   -DOUTPUT_ROOT=\"/your/output/parent/\"
#ifndef OUTPUT_ROOT
#define OUTPUT_ROOT "./results/"
#endif
#define BASE_PATH OUTPUT_ROOT

#define RUN_FOLDER "fig5_rerun_A0_0_2p25_5"
#define PATH BASE_PATH RUN_FOLDER "/"

static const int n_fig5_snapshot_times = 4;
static const double fig5_snapshot_times[n_fig5_snapshot_times] = {0.0, 10.0, 50.0, 100.0};
static const double fig5_A0_values[] = {0.0, 2.25, 5.0};
static const int n_fig5_A0_values = sizeof(fig5_A0_values) / sizeof(fig5_A0_values[0]);

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

    // Remove a trailing slash except for the root directory "/".
    if (len > 1 && tmp[len - 1] == '/') {
        tmp[len - 1] = '\0';
    }

    for (char* p = tmp + 1; *p; ++p) {
        if (*p == '/') {
            *p = '\0';
            ensure_dir(tmp);
            *p = '/';
        }
    }
    ensure_dir(tmp);
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
    double first_term = 0.0, second_term = 0.0, third_term = 0.0, fourth_term = 0.0, fifth_term = 0.0;

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

static void write_field_file(const char* field_file,
                             double **A, double **F, double **R) {
    FILE *fsi = fopen(field_file, "w");
    if (!fsi) {
        perror("fail to open field file");
        exit(1);
    }

    // x, y, A, F, R
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            fprintf(fsi, "%lf %lf %.12lf %.12lf %.12lf\n",
                    i * dx, j * dy, A[i][j], F[i][j], R[i][j]);
        }
        fprintf(fsi, "\n");
    }

    fclose(fsi);
}

static void save_field_snapshot_by_index(const char* dir_data, int z,
                                         double **A, double **F, double **R) {
    char field_file[1024];
    snprintf(field_file, sizeof(field_file), "%s/field_%04d.dat", dir_data, z);
    write_field_file(field_file, A, F, R);
}

static void save_field_snapshot_by_time(const char* dir_data, double t,
                                        double **A, double **F, double **R) {
    char field_file[1024];
    int t_int = static_cast<int>(llround(t));
    snprintf(field_file, sizeof(field_file), "%s/field_t%03d.dat", dir_data, t_int);
    write_field_file(field_file, A, F, R);
}

static void initialize_fig5_fungal_patch(double **F) {
    // Figure 5 uses the physical-domain initial condition described in the manuscript:
    // center (5,5), radius 2, F=0.5 inside, F=0 outside.
    // This deliberately avoids index-based radius mistakes such as "radius = 2 cells".
    const double cx = 5.0;
    const double cy = 5.0;
    const double radius = 2.0;
    const double inside = 0.5;
    const double outside = 0.0;
    const double r2 = radius * radius;

    for (int i = 0; i < Nx; i++) {
        const double x = i * dx;
        for (int j = 0; j < Ny; j++) {
            const double y = j * dy;
            const double dist2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            F[i][j] = (dist2 <= r2) ? inside : outside;
        }
    }
}


static double run_simulation(double A0, FILE* fSweep) {
    // -------------------------------------
    // Memory allocation
    // -------------------------------------
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

    // -------------------------------------
    // Diffusion coefficients for ADI
    // -------------------------------------
    double rAx = (dA * 0.5 * dt) / (dx * dx);
    double rAy = (dA * 0.5 * dt) / (dy * dy);
    double rFx = (dF * 0.5 * dt) / (dx * dx);
    double rFy = (dF * 0.5 * dt) / (dy * dy);

    // -------------------------------------
    // Initial condition
    // -------------------------------------
    srand(1);

    initialization(R, Nx, Ny, initR);
    initialization(A, Nx, Ny, A0);

    initialize_fig5_fungal_patch(F);

    // -------------------------------------
    // LU decomposition setup (Neumann BC)
    // -------------------------------------
    for (int i = 0; i < Nx; i++) {
        A_Ax[i] = 1 + 2 * rAx;
        B_Ax[i] = -rAx;
        C_Ax[i] = -rAx;
    }
    A_Ax[0] = 1 + 2 * rAx;
    B_Ax[0] = 0.0;
    C_Ax[0] = -2 * rAx;
    A_Ax[Nx - 1] = 1 + 2 * rAx;
    B_Ax[Nx - 1] = -2 * rAx;
    C_Ax[Nx - 1] = 0.0;

    for (int j = 0; j < Ny; j++) {
        A_Ay[j] = 1 + 2 * rAy;
        B_Ay[j] = -rAy;
        C_Ay[j] = -rAy;
    }
    A_Ay[0] = 1 + 2 * rAy;
    B_Ay[0] = 0.0;
    C_Ay[0] = -2 * rAy;
    A_Ay[Ny - 1] = 1 + 2 * rAy;
    B_Ay[Ny - 1] = -2 * rAy;
    C_Ay[Ny - 1] = 0.0;

    LUx(A_Ax, B_Ax, C_Ax, L_Ax, U_Ax, Nx);
    LUy(A_Ay, B_Ay, C_Ay, L_Ay, U_Ay, Ny);

    for (int i = 0; i < Nx; i++) {
        A_Fx[i] = 1 + 2 * rFx;
        B_Fx[i] = -rFx;
        C_Fx[i] = -rFx;
    }
    A_Fx[0] = 1 + 2 * rFx;
    B_Fx[0] = 0.0;
    C_Fx[0] = -2 * rFx;
    A_Fx[Nx - 1] = 1 + 2 * rFx;
    B_Fx[Nx - 1] = -2 * rFx;
    C_Fx[Nx - 1] = 0.0;

    for (int j = 0; j < Ny; j++) {
        A_Fy[j] = 1 + 2 * rFy;
        B_Fy[j] = -rFy;
        C_Fy[j] = -rFy;
    }
    A_Fy[0] = 1 + 2 * rFy;
    B_Fy[0] = 0.0;
    C_Fy[0] = -2 * rFy;
    A_Fy[Ny - 1] = 1 + 2 * rFy;
    B_Fy[Ny - 1] = -2 * rFy;
    C_Fy[Ny - 1] = 0.0;

    LUx(A_Fx, B_Fx, C_Fx, L_Fx, U_Fx, Nx);
    LUy(A_Fy, B_Fy, C_Fy, L_Fy, U_Fy, Ny);

    // -------------------------------------
    // Output folder preparation
    // -------------------------------------
    ensure_dir_recursive(PATH);
    ensure_dir_recursive(PATH "data");
    ensure_dir_recursive(PATH "region");

    char dir_data[1024];
    char dir_region[1024];

    snprintf(dir_data, sizeof(dir_data), PATH "data/A0_%.3f", A0);
    snprintf(dir_region, sizeof(dir_region), PATH "region/A0_%.3f", A0);

    ensure_dir_recursive(dir_data);
    ensure_dir_recursive(dir_region);

    char region_file[1024];
    snprintf(region_file, sizeof(region_file), "%s/region.txt", dir_region);

    FILE *fSum = fopen(region_file, "w");
    if (!fSum) {
        perror("fail to open region file");
        exit(1);
    }

    char snapshot_index_file[1024];
    snprintf(snapshot_index_file, sizeof(snapshot_index_file), "%s/snapshot_index.txt", dir_data);
    FILE *fSnap = fopen(snapshot_index_file, "w");
    if (!fSnap) {
        perror("fail to open snapshot index file");
        exit(1);
    }

    fprintf(fSum, "# time total_A total_F fungal_area_percent\n");
    fprintf(fSnap, "# time filename total_A total_F fungal_area_percent\n");

    tbb::task_arena arena;

    double saved_area_percent[n_eval_times] = {0.0};
    double saved_total_A[n_eval_times] = {0.0};
    double saved_total_F[n_eval_times] = {0.0};
    bool saved_once[n_eval_times] = {false};

    int save_step[n_eval_times];

    for (int q = 0; q < n_eval_times; q++) {
        save_step[q] = static_cast<int>(eval_times[q] / dt + 0.5);
    }

    int fig5_snapshot_step[n_fig5_snapshot_times];
    bool fig5_snapshot_saved[n_fig5_snapshot_times] = {false};
    for (int q = 0; q < n_fig5_snapshot_times; q++) {
        fig5_snapshot_step[q] = static_cast<int>(fig5_snapshot_times[q] / dt + 0.5);
    }

    // Save the exact initial condition before any time step is taken.
    // This is essential for Figure 5: all three A0 rows must have the same
    // large circular patch at t=0.
    double total_A0 = compute_total_mass(A);
    double total_F0 = compute_total_mass(F);
    double area_percent0 = fungal_area_percent(F);
    fprintf(fSum, "%.12lf %.12lf %.12lf %.12lf\n", 0.0, total_A0, total_F0, area_percent0);
    save_field_snapshot_by_time(dir_data, 0.0, A, F, R);
    fprintf(fSnap, "%.12lf field_t%03d.dat %.12lf %.12lf %.12lf\n",
            0.0, 0, total_A0, total_F0, area_percent0);
    fflush(fSum);
    fflush(fSnap);

    for (int q = 0; q < n_fig5_snapshot_times; q++) {
        if (fig5_snapshot_step[q] == 0) {
            fig5_snapshot_saved[q] = true;
        }
    }
    
    // -------------------------------------
    // Time iteration
    // -------------------------------------
    for (int k = 0; k < iter; k++)
    {
        // Neumann BC
        for (int j = 0; j < Ny; j++)
        {
            A[0][j]      = A[1][j];
            A[Nx - 1][j] = A[Nx - 2][j];
            F[0][j]      = F[1][j];
            F[Nx - 1][j] = F[Nx - 2][j];
        }

        for (int i = 0; i < Nx; i++)
        {
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

                        Fx[j][i] = coefficient_dF(F[i][j], rFy,
                                                  F[i][jp], F[i][j], F[i][jm],
                                                  R[i][j], m, theta);
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

                        Fy[i][j] = coefficient_dF(F_new[i][j], rFx,
                                                  F_new[ip][j], F_new[i][j], F_new[im][j],
                                                  R_new[i][j], m, theta);
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

        const int step = k + 1;
        const double current_time = step * dt;

        double total_A = compute_total_mass(A);
        double total_F = compute_total_mass(F);
        double area_percent = fungal_area_percent(F);

        for (int q = 0; q < n_eval_times; q++) {
            if (!saved_once[q] && step >= save_step[q]) {
                saved_total_A[q] = total_A;
                saved_total_F[q] = total_F;
                saved_area_percent[q] = area_percent;
                saved_once[q] = true;

                printf("Saved summary at t = %.6lf: total_F = %.12lf, area%% = %.12lf\n",
                    eval_times[q], saved_total_F[q], saved_area_percent[q]);
            }
        }

        for (int q = 0; q < n_fig5_snapshot_times; q++) {
            if (!fig5_snapshot_saved[q] && step >= fig5_snapshot_step[q]) {
                const double target_t = fig5_snapshot_times[q];
                const int target_t_int = static_cast<int>(llround(target_t));

                save_field_snapshot_by_time(dir_data, target_t, A, F, R);
                fprintf(fSnap, "%.12lf field_t%03d.dat %.12lf %.12lf %.12lf\n",
                        target_t, target_t_int, total_A, total_F, area_percent);
                fflush(fSnap);

                fig5_snapshot_saved[q] = true;
                printf("Saved Figure 5 snapshot: A0 = %.3lf, t = %.6lf, area%% = %.6lf\n",
                       A0, target_t, area_percent);
            }
        }

        if (step % cut == 0) {
            fprintf(fSum, "%.12lf %.12lf %.12lf %.12lf\n",
                    current_time, total_A, total_F, area_percent);
            fflush(fSum);

            printf("A0: %.3lf time: %.6lf area%%: %.6lf\n",
                   A0, current_time, area_percent);

            if (F[Nx / 2][Ny / 2] > 1.0 || isnan(F[Nx / 2][Ny / 2])) {
                printf("Unstable at step %d (t = %.6lf): F = %.6lf\n",
                       step, current_time, F[Nx / 2][Ny / 2]);
                break;
            }
        }
    }

    fclose(fSum);
    fclose(fSnap);

    // sweep summary
    for (int q = 0; q < n_eval_times; q++) {
        if (!saved_once[q]) {
            double total_A = compute_total_mass(A);
            double total_F = compute_total_mass(F);
            double area_percent = fungal_area_percent(F);

            saved_total_A[q] = total_A;
            saved_total_F[q] = total_F;
            saved_area_percent[q] = area_percent;

            printf("Warning: target time t = %.6lf not reached. Saved final value instead.\n", eval_times[q]);
        }
    }
    fprintf(fSweep, "%.12lf", A0);

    for (int q = 0; q < n_eval_times; q++) {
        fprintf(fSweep, " %.12lf %.12lf", saved_total_F[q], saved_area_percent[q]);
    }

    fprintf(fSweep, "\n");
    fflush(fSweep);

    // free memory
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

    return saved_area_percent[1];
}

int main(int argc, char **argv)
{
    printf("\n=====================================\n");
    printf("Figure 5 rerun batch\n");
    printf("Output root: %s\n", PATH);
    printf("=====================================\n");

    ensure_dir_recursive(PATH);
    ensure_dir_recursive(PATH "data");
    ensure_dir_recursive(PATH "region");
    ensure_dir_recursive(PATH "summary");

    char summary_file[1024];
    snprintf(summary_file, sizeof(summary_file), PATH "summary/summary_fig5_rerun.txt");

    FILE *fSweep = fopen(summary_file, "w");
    if (!fSweep) {
        perror("fail to open summary file");
        return 1;
    }

    fprintf(fSweep, "# A0");
    for (int q = 0; q < n_eval_times; q++) {
        fprintf(fSweep, " total_F_t%.6g area_t%.6g", eval_times[q], eval_times[q]);
    }
    fprintf(fSweep, "\n");
    fflush(fSweep);

    std::vector<double> doses;
    if (argc > 1) {
        // Optional override: ./program 0 2.25 5
        for (int i = 1; i < argc; i++) {
            doses.push_back(atof(argv[i]));
        }
    } else {
        for (int i = 0; i < n_fig5_A0_values; i++) {
            doses.push_back(fig5_A0_values[i]);
        }
    }

    for (double A0 : doses) {
        printf("\n=====================================\n");
        printf("Running simulation : A0 = %.6lf\n", A0);
        printf("=====================================\n");

        run_simulation(A0, fSweep);

        printf("Simulation finished for A0 = %.6lf\n", A0);
    }

    fclose(fSweep);

    printf("\nAll Figure 5 rerun simulations finished.\n");
    printf("Saved snapshots are in:\n");
    printf("  %sdata/A0_0.000/field_t000.dat, field_t010.dat, field_t050.dat, field_t100.dat\n", PATH);
    printf("  %sdata/A0_2.250/field_t000.dat, field_t010.dat, field_t050.dat, field_t100.dat\n", PATH);
    printf("  %sdata/A0_5.000/field_t000.dat, field_t010.dat, field_t050.dat, field_t100.dat\n", PATH);

    return 0;
}
