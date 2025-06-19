# README - MALDI-GeoMx-Spatial-Transcriptomics Nitrogen Index calculation with IsotopeakeR

This repository implements two complementary spectral processing pipelines for MALDI imaging data:

- **Nitrogen Index RP** in reflector mode: processes raw spectra, applies SG smoothing, SNIP baseline correction, MAD- and valley-based peak detection, computes top-N peak metrics and ratios, then exports tables and overview plots.
- **Nitrogen Index LP** in linear positive mode: extracts a targeted m/z region (amyloid Aβ₁₋₄₂), performs normalization, smoothing, ALS baseline correction, FWHM peak analysis, and exports annotated results.

Use the sections below to explore each workflow in detail.

<details>
  
<summary><strong> Nitrogen Index RP (Reflector mode in MALDI MSI) </strong></summary>

---

### Table of Contents

1. [Step 1: Parameter Setup](#step-1-parameter-setup)  
2. [Step 2: Helper Functions](#step-2-helper-functions)  
3. [Step 3: Single-File Processing](#step-3-single-file-processing)  
4. [Step 4: Intensity Scaling](#step-4-intensity-scaling)  
5. [Step 5: Smoothing & Baseline Removal](#step-5-smoothing--baseline-removal)  
6. [Step 6: Peak Detection within ROI](#step-6-peak-detection-within-roi)  
7. [Step 7: Peak Metrics Computation](#step-7-peak-metrics-computation)  
8. [Step 8: Plot Generation](#step-8-plot-generation)  
9. [Step 9: Aggregation & Ratio Calculation](#step-9-aggregation--ratio-calculation)  
10. [Step 10: Export & Visualization Output](#step-10-export--visualization-output)  

---

#### Step 1: Parameter Setup  
All tunable settings at top of **pipeline.R**:
- **folder_in**, **file_pat** – input directory & file pattern  
- **spec_min**, **spec_max** – spectrum ROI m/z bounds  
- **peak_min**, **peak_max** – peak-detection ROI bounds  
- **sgolay_p**, **sgolay_n** – SG filter order & window  
- **snip_iter** – SNIP baseline iterations  
- **SNR**, **halfWindow**, **minIntensity** – detection thresholds  
- **intensity_transform** – `"sqrt"`  
- **top_n** – how many peaks to keep  
- **scaling** – `"RMSE"` or `"pct"`  
- **ratio_from**, **ratio_to**, **ratio_basis** – peak numbers & basis  
- **out_table**, **out_plots** – Excel & PNG filenames  

---

#### Step 2: Helper Functions  
- **load_file(path)** – read TXT/CSV, numeric mz & intensity  
- **baseline_snip(y,it)** – iterative SNIP baseline  
- **detect_peaks_mad(...)** – local maxima + MAD filter  
- **detect_peaks_min_distance(...)** – min‐distance peak filter  
- **find_bounds_valley(x,y,idx)** – valley boundaries  
- **peak_metrics(x,y,L,R)** – AUC, centroid, FWHM  
- **transform_intensity(df,method)** – none/√/log transform  

---

#### Step 3: Single-File Processing  
**process_single(path)**:
1. Load & filter `[spec_min,spec_max]`  
2. Transform intensity  
3. Scale (RMSE or pct)  
4. SG smoothing + SNIP baseline removal  
5. Detect peaks in `[peak_min,peak_max]`, keep top_n  
6. Compute peak metrics & build tibble  
7. Return list(table, plot, spec)  

---

#### Step 4: Intensity Scaling  
- **RMSE**: divide by √mean(intensity²)  
- **Pct**: normalize max(intensity)=100  

---

#### Step 5: Smoothing & Baseline Removal  
1. **Savitzky–Golay** (sgolay_p, sgolay_n)  
2. **Subtract** SNIP baseline (snip_iter)  
3. **Clamp** negatives to zero  

---

#### Step 6: Peak Detection within ROI  
- Subset to `[peak_min,peak_max]`  
- `detect_peaks_min_distance()` with min_distance=0.5  
- Keep up to **top_n** peaks  

---

#### Step 7: Peak Metrics Computation  
For each peak:
- Boundaries via `find_bounds_valley()`  
- AUC/centroid/FWHM via `peak_metrics()`  
- Assemble file, peak#, mz, intensity, range, metrics  

---

#### Step 8: Plot Generation  
- `ggplot2` line plot of processed spectrum  
- Highlight ROIs, annotate peaks `#1`, `#2`, …  
- Title = filename  

---

#### Step 9: Aggregation & Ratio Calculation  
1. `list.files()` → `proc_list`  
2. `bind_rows()` → `metrics_all`  
3. `calc_peak_ratio_info()` → `ratio_table`  

---

#### Step 10: Export & Visualization Output  
- **Excel**: sheets “Peaks” & “Ratio” via `openxlsx` → `out_table`  
- **PNG grid**: arrange all plots via `gridExtra` → `out_plots`  
- Console message: n peaks in n spectra processed  

---
</details>

<details>
<summary><strong> Nitrogen Index LP (Linear positive mode in MALDI MSI) </strong></summary>

---

### Table of Contents

1. [Step 1: Creation of Raw Mass Spectra](#step-1-creation-of-raw-mass-spectra)  
2. [Step 2: Intensity Normalization](#step-2-intensity-normalization)  
3. [Step 3: Spectral Smoothing](#step-3-spectral-smoothing)  
4. [Step 4: Baseline Correction](#step-4-baseline-correction)  
5. [Step 5: Data Saving](#step-5-data-saving)  
6. [Step 6: Peak Analysis](#step-6-peak-analysis)  
7. [Step 7: Visualization](#step-7-visualization)  
8. [Step 8: Data Aggregation & Annotation](#step-8-data-aggregation--annotation)  
9. [Step 9: Sorting & Final Preparation](#step-9-sorting--final-preparation)  
10. [Step 10: Export of Results](#step-10-export-of-results)  

---

#### Step 1: Creation of Raw Mass Spectra

Raw mass spectra were generated by reading each input data file and extracting the signal around a targeted *m/z* value corresponding to amyloid Aβ₁₋₄₂. A ±25 *m/z* tolerance window was applied to isolate the region of interest before further processing.

---

#### Step 2: Intensity Normalization

Each spectrum’s intensity trace was scaled to percentage of its own maximum intensity. This “100% normalization” ensures comparability across spectra by compensating for absolute intensity differences between acquisitions.

---

#### Step 3: Spectral Smoothing  

To reduce noise while preserving peak shapes, a **Savitzky–Golay filter** was applied with:  
- Polynomial order 
- Window size  

This smoothing step enhances subsequent peak detection performance.

---

#### Step 4: Baseline Correction 

Smoothed spectra were corrected for low‐frequency background using the **Asymmetric Least Squares (ALS)** algorithm with:  
- Smoothing parameter (λ)
- Asymmetry parameter (p)  
- Maximum iterations  

Baseline drift is thereby removed, improving accuracy of peak quantification.

---

#### Step 5: Data Saving

Baseline-corrected spectra are serialized to `.rds` files. This allows you to resume downstream analyses without repeating the computationally intensive preprocessing steps.

---

#### Step 6: Peak Analysis  
Automatic peak analysis is performed on each baseline-corrected spectrum to extract quantitative descriptors. In this pipeline we focus on:  
- **Full Width at Half Maximum (FWHM)** — measures peak broadness.  

Additional metrics can be added as needed.

---

#### Step 7: Visualization  

Analysis results are plotted using a grid layout that arranges individual spectrum overlays side by side. This layout facilitates rapid visual comparison of peak shapes and heights across all samples.

---

#### Step 8: Data Aggregation & Annotation  

Individual peak tables are concatenated into a single data frame. Row names or file identifiers are parsed to extract metadata fields—specifically **model type** and **region of interest (ROI)**—which are added as separate columns for structured annotation.

---

#### Step 9: Sorting & Final Preparation  

The combined dataset is sorted by numeric ROI identifier. This ordering ensures that results are organized in biologically or spatially meaningful sequence prior to export.

---

#### Step 10: Export of Results  

The finalized, annotated dataset was exported to an **Excel spreadsheet** to facilitate downstream statistical analysis, reporting, and sharing with collaborators.

---
</details>


