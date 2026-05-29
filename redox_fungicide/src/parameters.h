#pragma once

// ============================================================
//  Model and discretization parameters
// ============================================================

// ===== Model parameters =====

// diffusion coefficients
inline constexpr double dA = 0.005;   // d_A = D_A / (sigma0 * L^2)
inline constexpr double dF = 0.001;   // d_F = D_F / (sigma0 * L^2)

// interaction parameters
inline constexpr double alpha = 50.0; // alpha = kappa * F_max / sigma0
inline constexpr double l     = 0.3;  // l = delta / sigma0
inline constexpr double eps   = 0.05; // epsilon = eta / sigma0
inline constexpr double beta  = 10.0; // beta = omega0 / sigma0

// nonlinear response parameters
inline constexpr double m      = 2.0;  // Hill exponent
inline constexpr double theta  = 0.30; // threshold parameter


// ============================================================
//  >>> CONVERGENCE-STUDY BLOCK (Tier A item A1) <<<
//
//  For the grid / time-step convergence study, change ONLY the
//  values in this block and rebuild (./exec.sh). Nothing else
//  needs to be edited: dx, dy are derived from Nx, Ny, and the
//  time step dt is picked up automatically by functions.cpp.
//
//  SPATIAL refinement : change Nx, Ny  (keep dt fixed)
//      e.g.  Nx = Ny = 125, 250, 500
//  TEMPORAL refinement: change dt      (keep Nx, Ny fixed)
//      and set iter = (t_final / dt) + 1   so t_final stays 200
//      e.g.  dt = 4e-4 -> iter =  500001
//            dt = 2e-4 -> iter = 1000001
//            dt = 1e-4 -> iter = 2000001
// ============================================================

// ----- Spatial domain -----
inline constexpr double Lx = 10.0;
inline constexpr double Ly = 10.0;

inline constexpr int Nx = 500;     // <-- change for spatial convergence
inline constexpr int Ny = 500;     // <-- change for spatial convergence

inline constexpr double dx = Lx / (Nx - 1);   // derived; do not edit
inline constexpr double dy = Ly / (Ny - 1);   // derived; do not edit

// ----- Time stepping -----
inline constexpr double dt   = 0.0001;        // <-- change for temporal convergence
inline constexpr int    iter = 2000001;       // <-- set to (t_final / dt) + 1 ; t_final = 200
inline constexpr int    cut  = 1000;          // output stride (in steps)

// ----- Evaluation times (used for the summary and convergence) -----
inline constexpr int    n_eval_times = 4;
inline constexpr double eval_times[n_eval_times] = {20.0, 50.0, 100.0, 200.0};

// ============================================================
//  >>> END OF CONVERGENCE-STUDY BLOCK <<<
// ============================================================


// ===== Initial conditions =====
inline constexpr double initR = 1.0;

// ===== Fungal infection patch (medium patch) =====
inline constexpr double F_center_cx = 5.0;
inline constexpr double F_center_cy = 5.0;
inline constexpr double F_radius    = 2.0;    // medium-sized patch
inline constexpr double F_inside    = 0.5;   // starts below F_infty
inline constexpr double F_outside   = 0.0;

inline constexpr double F_area_threshold = 0.001;

// ===== Noise =====
inline constexpr double A_noise = 0.01;
inline constexpr double F_noise = 0.01;

// ===== Pattern analysis =====
inline constexpr int Lp = 2;
