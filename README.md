# README - MALDI-GeoMx-Spatial-Transcriptomics Nitrogen Index Calculation with IsotopeakeR

This repository implements a web-based Shiny application, **IsotopeakeR Beta**, for interactive preprocessing and isotopic analysis of MALDI mass spectrometry imaging (MSI) data, alongside a spectral processing pipeline for linear positive mode:

- **IsotopeakeR Web App for Nitrogen Index RP**: A browser-based tool that allows users to upload raw spectra (.dat, .txt, .csv), apply preprocessing steps (subsetting, intensity transformation, Savitzky-Golay smoothing, SNIP baseline correction, TIC/RMSE normalization), detect peaks using MAD- and valley-based methods, visualize results with interactive `plotly` plots, and export peak data as CSV.
- **Nitrogen Index LP** (Linear Positive Mode): Extracts a targeted m/z region (e.g., amyloid Aβ₁₋₄₂), performs normalization, smoothing, ALS baseline correction, FWHM peak analysis, and exports annotated results.

Use the sections below to explore the web application and the Nitrogen Index LP pipeline in detail.

## Overview

<details>
  
<summary><strong> IsotopeakeR Beta for Nitrogen Index RP </strong></summary>

---

**IsotopeakeR Beta** is a Shiny web application designed for preprocessing and isotopic analysis of MALDI MSI data. It supports `.dat`, `.txt`, and `.csv` file formats, offering tools for spectral preprocessing, peak detection, interactive visualization, and data export. Deployed on a web server (e.g., Shiny Server or shinyapps.io), it provides a user-friendly interface accessible via a browser, eliminating the need for local R installation.

**Author**: Maciej Dulewicz (maciej.dulewicz@gu.se)

## Features

- **Preprocessing**: Subset m/z ranges, apply sqrt/log intensity transformation, Savitzky-Golay smoothing, SNIP baseline correction, and TIC/RMSE normalization.
- **Peak Detection**: Identify peaks within a specified m/z range using MAD-based filtering with customizable SNR, window size, and minimum distance.
- **Visualization**: Interactive `plotly` plots for raw and processed spectra, with annotated peaks and integration boundaries.
- **Export**: Download peak data as CSV for further analysis.
- **Web Access**: Hosted on a server for browser-based usage, supporting file uploads up to 100 MB.

## Installation and Deployment

1. Access the app via the provided URL (e.g., https://maciejdulewiczgu.shinyapps.io/IsotopeakeR_Beta/).

## Usage Guide: IsotopeakeR Web Application

### Table of Contents
1. [Accessing the Web Application](#accessing-the-web-application)
2. [Preprocessing Panel Overview](#preprocessing-panel-overview)
3. [Uploading Spectra Files](#uploading-spectra-files)
4. [Configuring Preprocessing Options](#configuring-preprocessing-options)
5. [Processing Spectra](#processing-spectra)
6. [PeaksCatcher Panel Overview](#peakscatcher-panel-overview)
7. [Peak Detection Settings](#peak-detection-settings)
8. [Detecting and Analyzing Peaks](#detecting-and-analyzing-peaks)
9. [Visualizing Results](#visualizing-results)
10. [Exporting Peak Data](#exporting-peak-data)

#### Accessing the Web Application
- Open the app URL (e.g., https://maciejdulewiczgu.shinyapps.io/IsotopeakeR_Beta/) in a browser (Chrome, Firefox, Edge).
- No local R installation required; ensure a stable internet connection.
- The interface loads with a sandstone theme, titled "IsotopeakeR Beta – Preprocessing & Isotopic Analysis."

#### Preprocessing Panel Overview
- **Purpose**: Upload, preprocess, and visualize raw mass spectrometry data from https://maciejdulewiczgu.shinyapps.io/MALDI_GEOMX_VOLCANO/ or from Zenodo.
- **Components**:
  - **Sidebar**: File upload controls, preprocessing settings, and action buttons.
  - **Main Panel**: Interactive `plotly` spectrum preview and processing log.
- **Navigation**: Switch between spectra using "Previous," "Next," or the "Select RAW Spectrum" dropdown.

#### Uploading Spectra Files
- **Formats**: `.dat`, `.txt` (space-separated), or `.csv` with m/z and intensity columns.
- **Steps**:
  1. Click "Select raw spectra files" in the sidebar.
  2. Upload files (total size ≤ 100 MB).
  3. Click "Load files"; a pop-up confirms the number of files loaded.
  4. Clear uploaded spectra with "Clear spectra."
- **Result**: Spectra names appear in the dropdown for selection.

#### Configuring Preprocessing Options
- **Subset m/z Range**:
  - Enable "Process m/z sub-range?" to restrict processing.
  - Set `m/z min` (default: 3500) and `m/z max` (default: 5500).
- **Intensity Transformation**:
  - Check "Transform Intensity?"; select `sqrt` (default), `log`, or `none`.
- **Smoothing (Savitzky-Golay)**:
  - Enable "Smooth (Savitzky-Golay)?".
  - Set `Polynomial order` (default: 3) and `Window size` (default: 11, must exceed polynomial order).
- **Baseline Correction (SNIP)**:
  - Enable "Baseline Correction (SNIP)?".
  - Set `SNIP Iterations` (default: 100).
- **Normalization**:
  - Select `none` (default), `TIC` (Total Ion Current), or `RMSE` (Root Mean Square Error).

#### Processing Spectra
- **Single Spectrum**:
  1. Select a spectrum from the dropdown.
  2. Configure preprocessing options.
  3. Click "Process & Save (One)"; a pop-up shows the processed spectrum’s name (e.g., `original_proc_HHMMSS`).
- **All Spectra**:
  1. Set preprocessing options.
  2. Click "Process ALL & Save"; a pop-up confirms the number of processed spectra.
- **Reset**: Clear processed data with "Reset Processed Results."
- **Preview**: The `plotly` plot updates dynamically to reflect preprocessing settings.

#### PeaksCatcher Panel Overview
- **Purpose**: Detect and analyze peaks in processed spectra within a specified m/z range.
- **Components**:
  - **Sidebar**: Spectra selection, peak detection parameters, and export button.
  - **Main Panel**: Interactive `plotly` plot, peak summary table, spectra details, and statistical summary.
- **Access**: Navigate to the "PeaksCatcher" tab after processing spectra.

#### Peak Detection Settings
- **Spectra Selection**: Select one or more processed spectra (multiple selections enabled).
- **Search Range**:
  - Set `m/z Start` (default: 4500) and `m/z End` (default: 4600).
  - Optionally set `Manual Baseline` (default: 0) to subtract a constant intensity.
- **Peak Detection Parameters**:
  - `SNR` (default: 2): Signal-to-noise ratio threshold.
  - `Half window size` (default: 20): Window for local maxima detection.
  - `Min Distance (m/z)` (default: 0.8): Minimum m/z separation between peaks.
  - `Max peaks` (default: 5): Maximum number of peaks to detect per spectrum.

#### Detecting and Analyzing Peaks
- **Action**: Click "Detect Peaks" to analyze selected spectra.
- **Process**:
  1. Subset spectra to the specified m/z range.
  2. Apply manual baseline subtraction.
  3. Detect peaks using MAD-based filtering (`my_detect_peaks()`).
  4. Identify peak boundaries with valley-based method (`find_bounds_valley()`).
  5. Compute metrics (AUC, centroid, FWHM) using `calc_auc()` and `calc_peak_params()`.
- **Output**:
  - **Peak Summary Table**: Peak details (spectrum, m/z, intensity, AUC, centroid, FWHM).
  - **Spectra Summary Table**: Number of peaks detected per spectrum.
  - **Stats Table**: Average and median m/z per spectrum.

#### Visualizing Results
- **Plot**:
  - Interactive `plotly` plot displays spectra in the m/z range.
  - Peaks marked with red dots, labeled with peak number and m/z.
  - Green dotted lines show peak integration boundaries.
  - Search range shaded in light blue.
- **Interactivity**: Hover for peak details; zoom/pan with `plotly` controls.

#### Exporting Peak Data
- **Action**: Click "Export Peak Data (CSV)" to download peak summary.
- **Output**: `.csv` file (e.g., `peak_data_YYYY-MM-DD.csv`) with spectrum name, peak number, m/z, intensity, range, AUC, centroid, and FWHM.
- **Use Case**: Import into Excel or R for further analysis.

---

</details>

<details>
<summary><strong>Nitrogen Index LP (Linear Positive Mode in MALDI MSI)</strong></summary>

### Table of Contents
1. [Creation of Raw Mass Spectra](#step-1-creation-of-raw-mass-spectra)
2. [Intensity Normalization](#step-2-intensity-normalization)
3. [Spectral Smoothing](#step-3-spectral-smoothing)
4. [Baseline Correction](#step-4-baseline-correction)
5. [Data Saving](#step-5-data-saving)
6. [Peak Analysis](#step-6-peak-analysis)
7. [Visualization](#step-7-visualization)
8. [Data Aggregation & Annotation](#step-8-data-aggregation--annotation)
9. [Sorting & Final Preparation](#step-9-sorting--final-preparation)
10. [Export of Results](#step-10-export-of-results)

#### Step 1: Creation of Raw Mass Spectra
Raw mass spectra are generated by reading input files and extracting the signal around a targeted m/z value (e.g., amyloid Aβ₁₋₄₂) with a ±25 m/z tolerance window.

#### Step 2: Intensity Normalization
Spectra are scaled to percentage of maximum intensity (100% normalization) to ensure comparability across acquisitions.

#### Step 3: Spectral Smoothing
Apply **Savitzky-Golay filter** with:
- Polynomial order.
- Window size.
This reduces noise while preserving peak shapes.

#### Step 4: Baseline Correction
Use **Asymmetric Least Squares (ALS)** with:
- Smoothing parameter (λ).
- Asymmetry parameter (p).
- Maximum iterations.
Removes baseline drift for accurate peak quantification.

#### Step 5: Data Saving
Serialize baseline-corrected spectra to `.rds` files for efficient downstream analysis.

#### Step 6: Peak Analysis
Perform automatic peak analysis to extract:
- **Full Width at Half Maximum (FWHM)**: Measures peak broadness.
Additional metrics can be included as needed.

#### Step 7: Visualization
Plot results in a grid layout for side-by-side comparison of peak shapes and heights across samples.

#### Step 8: Data Aggregation & Annotation
Concatenate peak tables into a single data frame, parsing file identifiers for metadata (e.g., model type, ROI).

#### Step 9: Sorting & Final Preparation
Sort the dataset by numeric ROI identifier for biologically/spatially meaningful organization.

#### Step 10: Export of Results
Export the annotated dataset to an **Excel spreadsheet** for downstream analysis and sharing.

</details>

## Issues
Report bugs or suggest features via the [GitHub Issues](https://github.com/yourusername/isotopeaker/issues) page.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact
For questions, contact Maciej Dulewicz at maciej.dulewicz@gu.se.
