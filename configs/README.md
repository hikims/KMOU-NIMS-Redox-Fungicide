# Configuration files

| File | Used by | Paper item |
|---|---|---|
| `b2_cases.txt` | `main_fig2_table3_4_dose_time_maps` | Main Figure 2 and Table 2 (fixed-time critical doses) |
| `threshold_list.txt` | `main_supp_figS1_threshold` | Multiplicative factors on dA, dF, alpha for the parameter-dependence scan |
| `list.txt` | `main_supp_figS1_threshold` | Dose scan feeding Supplementary Figures S1, S2 and S3 |

These files match the manuscript settings. Dose ranges and grid resolution may
be adjusted before rerunning the full simulations.

## Note on the dose range in `list.txt`

`list.txt` is the scan actually used for the submitted Supplementary Figures
S1-S3. It was regenerated from
`results/threshold_parameter/summary/threshold_scan_results.txt`, so rerunning
the driver reproduces those figures.

`A0` starts at 1.0 and advances in steps of 0.1, but the upper end differs per
parameter setting because each scan stopped once the suppression criterion was
met: 3.1 for `dF` x 0.5, 3.2 for `dA` x 0.5, 3.5 for `dF` x 0.25, 3.6 for
`alpha` x 1.5, 4.0 for `dA` x 0.25, 4.8 for `alpha` x 2, and 6.0 for the
remaining eight settings. There are 559 cases in total.

## Known issue in the output header

The driver writes a 16-name header to `threshold_scan_results.txt`,

```
# tag group factor A0 dA_factor dF_factor alpha_factor dA dF alpha ...
```

but each data row carries 14 fields: `group` and `factor` are never written.
Any tool that maps the header names onto the columns is therefore misaligned by
two positions. `plot_figS1.py` and `plot_figS2_S3.py` ignore the header and
read by fixed position, so they are unaffected, but the header should be
corrected in the driver.
