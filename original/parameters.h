// #pragma once

// // ===== Model parameters =====

// // diffusion coefficients
// inline constexpr double dA = 0.005;   // d_A = D_A / (sigma0 * L^2)
// inline constexpr double dF = 0.001;   // d_F = D_F / (sigma0 * L^2)

// // interaction parameters
// inline constexpr double alpha = 50.0; // alpha = kappa * F_max / sigma0
// inline constexpr double l     = 0.3; // l = delta / sigma0
// inline constexpr double eps   = 0.05; // epsilon = eta / sigma0
// inline constexpr double beta  = 10.0; // beta = omega0 / sigma0

// // nonlinear response parameters
// inline constexpr double m      = 2.0; // Hill exponent
// inline constexpr double theta  = 0.30; // threshold parameter


// // ===== Spatial domain =====
// inline constexpr double Lx = 10.0;
// inline constexpr double Ly = 10.0;

// inline constexpr int Nx = 500;
// inline constexpr int Ny = 500;

// inline constexpr double dx = Lx / (Nx - 1);
// inline constexpr double dy = Ly / (Ny - 1);


// // ===== Time stepping =====
// inline constexpr double dt = 0.0001;
// inline constexpr int iter  = 1000001;
// inline constexpr int cut   = 1000;
// inline constexpr int n_eval_times = 3;
// inline constexpr double eval_times[n_eval_times] = {20.0, 50.0, 100.0};


// // ===== Initial conditions =====
// // inline constexpr double initA = 0.5;
// inline constexpr double initR = 1.0;

// // ===== Fungal infection patch =====
// inline constexpr double F_center_cx = 5.0;
// inline constexpr double F_center_cy = 5.0;
// inline constexpr double F_radius = 0.25;
// inline constexpr double F_inside = 0.5;
// inline constexpr double F_outside = 0.0;

// inline constexpr double F_area_threshold = 0.001;

// // ===== Noise =====
// inline constexpr double A_noise = 0.01;
// inline constexpr double F_noise = 0.01;

// // ===== Pattern analysis =====
// inline constexpr int Lp = 2;

#pragma once

// ===== Model parameters =====

// diffusion coefficients
inline constexpr double dA = 0.005;   // d_A = D_A / (sigma0 * L^2)
inline constexpr double dF = 0.001;   // d_F = D_F / (sigma0 * L^2)

// interaction parameters
inline constexpr double alpha = 50.0; // alpha = kappa * F_max / sigma0
inline constexpr double l     = 0.3;  // l = delta / sigma0
inline constexpr double eps   = 0.05; // epsilon = eta / sigma0
inline constexpr double beta  = 10.0; // beta = omega0 / sigma0  (default; can be overridden via CLI)

// nonlinear response parameters
inline constexpr double m      = 2.0;  // Hill exponent
inline constexpr double theta  = 0.30; // threshold parameter


// ===== Spatial domain =====
inline constexpr double Lx = 10.0;
inline constexpr double Ly = 10.0;

inline constexpr int Nx = 500;
inline constexpr int Ny = 500;

inline constexpr double dx = Lx / (Nx - 1);
inline constexpr double dy = Ly / (Ny - 1);


// ===== Time stepping =====
inline constexpr double dt = 0.0001;
inline constexpr int iter  = 2000001;     // t_final = 200
inline constexpr int cut   = 1000;
inline constexpr int n_eval_times = 4;
inline constexpr double eval_times[n_eval_times] = {20.0, 50.0, 100.0, 200.0};


// ===== Initial conditions =====
inline constexpr double initR = 1.0;

// ===== Fungal infection patch (Option B: medium patch) =====
inline constexpr double F_center_cx = 5.0;
inline constexpr double F_center_cy = 5.0;
inline constexpr double F_radius    = 2.0;    // medium-sized patch (was 0.25)
inline constexpr double F_inside    = 0.5;   
inline constexpr double F_outside   = 0.0;

inline constexpr double F_area_threshold = 0.001;

// ===== Noise =====
inline constexpr double A_noise = 0.01;
inline constexpr double F_noise = 0.01;

// ===== Pattern analysis =====
inline constexpr int Lp = 2;