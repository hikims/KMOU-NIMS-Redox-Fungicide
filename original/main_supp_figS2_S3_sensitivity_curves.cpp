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

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#define PATH "/Users/mkim/Desktop/Research/fungicide_jh/"

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

static void save_field_snapshot(const char* dir_data, int z,
                                double **A, double **F, double **R) {
    char field_file[1024];
    snprintf(field_file, sizeof(field_file), "%s/field_%04d.dat", dir_data, z);

    FILE *fsi = fopen(field_file, "w");
    if (!fsi) {
        perror("fail to open field file");
        exit(1);
    }

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            fprintf(fsi, "%lf %lf %.12lf %.12lf %.12lf\n",
                    i * dx, j * dy, A[i][j], F[i][j], R[i][j]);
        }
        fprintf(fsi, "\n");
    }

    fclose(fsi);
}

static double run_simulation(
    const std::string& tag,
    double A0,
    double dA_factor,
    double dF_factor,
    double alpha_factor,
    double dA_run,
    double dF_run,
    double alpha_run,
    FILE* fSweep
) {
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
             F_radius,
             F_inside, F_outside);

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

    ensure_dir(PATH "sensitivity");
    ensure_dir(PATH "sensitivity/data");
    ensure_dir(PATH "sensitivity/region");

    char dir_data[1024];
    char dir_region[1024];
    char dir_movie[1024];

    snprintf(dir_data, sizeof(dir_data), PATH "sensitivity/data/%s", tag.c_str());
    snprintf(dir_region, sizeof(dir_region), PATH "sensitivity/region/%s", tag.c_str());
    snprintf(dir_movie, sizeof(dir_movie), PATH "sensitivity/data/%s/png", tag.c_str());

    ensure_dir(dir_data);
    ensure_dir(dir_region);
    ensure_dir(dir_movie);

    char region_file[1024];
    snprintf(region_file, sizeof(region_file), "%s/region.txt", dir_region);

    FILE *fSum = fopen(region_file, "w");
    if (!fSum) {
        perror("fail to open region file");
        exit(1);
    }

    fprintf(fSum, "# time total_A total_F fungal_area_percent\n");

    tbb::task_arena arena;

    double saved_area_percent[n_eval_times] = {0.0};
    double saved_total_A[n_eval_times] = {0.0};
    double saved_total_F[n_eval_times] = {0.0};
    bool saved_once[n_eval_times] = {false};

    int save_step[n_eval_times];

    for (int q = 0; q < n_eval_times; q++) {
        save_step[q] = static_cast<int>(eval_times[q] / dt + 0.5);
    }

    for (int k = 0; k < iter; k++) {
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

        double total_A = compute_total_mass(A);
        double total_F = compute_total_mass(F);
        double area_percent = fungal_area_percent(F);

        for (int q = 0; q < n_eval_times; q++) {
            if (!saved_once[q] && k >= save_step[q]) {
                saved_total_A[q] = total_A;
                saved_total_F[q] = total_F;
                saved_area_percent[q] = area_percent;
                saved_once[q] = true;

                printf("[%s] Saved summary at t = %.6lf: total_F = %.12lf, area%% = %.12lf\n",
                       tag.c_str(), eval_times[q], saved_total_F[q], saved_area_percent[q]);
            }
        }

        if (k % cut == 0) {
            fprintf(fSum, "%.12lf %.12lf %.12lf %.12lf\n",
                    k * dt, total_A, total_F, area_percent);
            fflush(fSum);

            printf("[%s] time: %.6lf file: %d area%%: %.6lf\n",
                   tag.c_str(), dt * k, k / cut, area_percent);

            int z = k / cut;
            save_field_snapshot(dir_data, z, A, F, R);

            if (F[Nx / 2][Ny / 2] > 1.0 || isnan(F[Nx / 2][Ny / 2])) {
                printf("Unstable at step %d (t = %.6lf): F = %.6lf\n",
                       k, k * dt, F[Nx / 2][Ny / 2]);
                break;
            }
        }
    }

    fclose(fSum);

    for (int q = 0; q < n_eval_times; q++) {
        if (!saved_once[q]) {
            double total_A = compute_total_mass(A);
            double total_F = compute_total_mass(F);
            double area_percent = fungal_area_percent(F);

            saved_total_A[q] = total_A;
            saved_total_F[q] = total_F;
            saved_area_percent[q] = area_percent;

            printf("Warning: target time t = %.6lf not reached. Saved final value instead.\n",
                   eval_times[q]);
        }
    }

    fprintf(fSweep, "%s %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf",
            tag.c_str(),
            A0,
            dA_factor,
            dF_factor,
            alpha_factor,
            dA_run,
            dF_run,
            alpha_run);

    for (int q = 0; q < n_eval_times; q++) {
        fprintf(fSweep, " %.12lf %.12lf", saved_total_F[q], saved_area_percent[q]);
    }

    fprintf(fSweep, "\n");
    fflush(fSweep);

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

    return saved_area_percent[n_eval_times - 1];
}

int main()
{
    printf("\n=====================================\n");
    printf("Running sensitivity analysis using list.txt\n");
    printf("=====================================\n");

    ensure_dir(PATH "sensitivity");
    ensure_dir(PATH "sensitivity/data");
    ensure_dir(PATH "sensitivity/region");
    ensure_dir(PATH "sensitivity/summary");

    char summary_file[1024];
    snprintf(summary_file, sizeof(summary_file),
             PATH "sensitivity/summary/sensitivity_results.txt");

    FILE *fSweep = fopen(summary_file, "w");
    if (!fSweep) {
        perror("fail to open sensitivity summary file");
        return 1;
    }

    fprintf(fSweep,
        "# tag A0 dA_factor dF_factor alpha_factor dA dF alpha "
        "total_F_t20 area_t20 total_F_t50 area_t50 total_F_t100 area_t100\n");

    std::ifstream fin(PATH "./code/list.txt");
    if (!fin.is_open()) {
        std::cerr << "Cannot open list.txt: " << PATH << "list.txt" << std::endl;
        fclose(fSweep);
        return 1;
    }

    std::string line;
    int run_count = 0;

    while (std::getline(fin, line)) {
        if (line.size() == 0) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);

        std::string tag;
        double A0;
        double dA_factor;
        double dF_factor;
        double alpha_factor;

        if (!(iss >> tag >> A0 >> dA_factor >> dF_factor >> alpha_factor)) {
            std::cerr << "Skipping invalid line: " << line << std::endl;
            continue;
        }

        double dA_run = dA * dA_factor;
        double dF_run = dF * dF_factor;
        double alpha_run = alpha * alpha_factor;

        printf("\n-------------------------------------\n");
        printf("Run %d: %s\n", run_count + 1, tag.c_str());
        printf("A0          = %.6lf\n", A0);
        printf("dA_factor   = %.6lf, dA = %.12lf\n", dA_factor, dA_run);
        printf("dF_factor   = %.6lf, dF = %.12lf\n", dF_factor, dF_run);
        printf("alpha_factor= %.6lf, alpha = %.12lf\n", alpha_factor, alpha_run);
        printf("-------------------------------------\n");

        run_simulation(
            tag,
            A0,
            dA_factor,
            dF_factor,
            alpha_factor,
            dA_run,
            dF_run,
            alpha_run,
            fSweep
        );

        run_count++;
    }

    fin.close();
    fclose(fSweep);

    printf("\n=====================================\n");
    printf("Sensitivity analysis finished.\n");
    printf("Total runs: %d\n", run_count);
    printf("Summary saved to:\n%s\n", summary_file);
    printf("=====================================\n");

    return 0;
}