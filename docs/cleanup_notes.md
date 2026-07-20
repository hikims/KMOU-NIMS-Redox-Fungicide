# Cleanup notes

1. Driver filenames are preserved exactly as provided.
2. Absolute local paths were replaced with a compile-time configurable `OUTPUT_ROOT` macro, defaulting to `./results/`.
3. The uploaded `main_supp_figS1_threshold.cpp` was identical to `main_fig5_spatial_dynamics.cpp`. For a functional reproducibility package, it was replaced under the same filename with the threshold-scan driver corresponding to Supplementary Figure S1.
4. No large raw simulation outputs are included.
5. Shared solver files are expected in `src/` before compilation.

## Second pass (Figure 6 and path portability)

6. Added the Figure 6 pipeline: `drivers/main_fig6_front_speed.cpp`,
   `scripts/front_speed_exec.sh`, `scripts/plot_fig6.py`. The driver uses the
   same `OUTPUT_ROOT` compile-time macro as the other drivers, and the shell
   script delegates the build to the Makefile.
7. Removed the absolute path `/Users/mkim/Desktop/Research/fungicide_jh` from
   five plotting scripts. All paths are now resolved relative to the
   repository root; `REDOX_FUNGICIDE_RESULTS` overrides the results tree.
8. `Makefile`: added the `main_fig6_front_speed` target and implemented
   `TBB_PREFIX`, which the README documented but the Makefile did not read.
9. `plot_fig3.py` and `plot_fig4.py` now also write PDF, matching the
   `\includegraphics{fig3.pdf}` calls in the manuscript. Every figure script
   emits PDF, EPS and PNG with TrueType fonts embedded.
10. Supplementary Figures S1, S2 and S3 are all produced from
    `threshold_parameter/summary/threshold_scan_results.txt`, written by the
    threshold driver. `main_supp_figS2_S3_sensitivity_curves.cpp` is retained
    but its output is no longer read by any plotting script.
11. Fixed a grid-spacing inconsistency. `src/parameters.h` defines
    `dx = Lx/(Nx-1) = 10/499 = 0.0200401` on a node-based grid, but
    `plot_fig5.py` and `plot_fig6.py` used `h = 0.02`, and `plot_fig5.py`
    additionally placed the field at cell centres rather than at nodes. Both
    scripts now derive the spacing as `LX/(NX-1)` and use node coordinates
    `x = i*dx`. In `plot_fig6.py` this affects only the fallback path, since
    the front radius is computed from `dx*dy` in the data-file header when
    that header is present; where the fallback applies, the front radius and
    the fitted speed were understated by 0.2 %.
12. Rebuilt `configs/list.txt`. The previous version contained only the `dF`
    and `alpha` groups and scanned `A0` from 1.0 to 3.0, so it could not
    reproduce Supplementary Figure S1, which also plots a `dA` curve and whose
    largest critical dose lies above 4. The file now covers all three groups
    over `A0 = 1.0 ... 5.0`. `scripts/make_list_from_results.py` regenerates
    the file directly from a completed scan when an exact match is wanted.
12. `configs/list.txt` was regenerated from the scan file that produced the
    submitted Supplementary Figures S1-S3. The previous version stopped at
    `A0 = 3.0` and could not reproduce them; the scan actually used reaches
    `A0 = 6.0`, with a per-setting upper limit. Duplicate `A0 = 3.00` rows,
    left by merging the original scan with its extension, were removed.
13. Noted, but not changed: the header of `threshold_scan_results.txt` names
    16 columns while the rows carry 14. The plotting scripts read by position
    and are unaffected.
