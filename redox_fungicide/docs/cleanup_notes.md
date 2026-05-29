# Cleanup notes

1. Driver filenames are preserved exactly as provided.
2. Absolute local paths were replaced with a compile-time configurable `OUTPUT_ROOT` macro, defaulting to `./results/`.
3. The uploaded `main_supp_figS1_threshold.cpp` was identical to `main_fig5_spatial_dynamics.cpp`. For a functional reproducibility package, it was replaced under the same filename with the threshold-scan driver corresponding to Supplementary Figure S1.
4. No large raw simulation outputs are included.
5. Shared solver files are expected in `src/` before compilation.
