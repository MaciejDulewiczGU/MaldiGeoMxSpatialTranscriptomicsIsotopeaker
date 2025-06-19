#Instructions
#1. Upload .dat or .csv files containing m/z and intensity data using the 'Select raw spectra files' input.

#2. Click 'Load files' to import the spectra. Use 'Previous' and 'Next' to navigate between loaded spectra.

#3. Select preprocessing options (subset range, intensity transformation, smoothing, baseline correction, normalization).

#4. Click 'Process & Save (One)' to process the selected spectrum or 'Process ALL & Save' to process all loaded spectra.

#5. Use 'Reset Processed Results' to clear all processed data.


# ------------------------------------------------------------------------------
# FUNCTION DESCRIPTIONS – Helper Functions for Spectrum Processing Pipeline
# ------------------------------------------------------------------------------

# create_spectrum(mz_vec, intensity_vec, file_name='')
# Purpose:        Creates a spectrum object as a list.
# Inputs:         mz_vec – Numeric vector of m/z values.
#                 intensity_vec – Numeric vector of intensity values.
#                 file_name – String, name of the source file (default: empty).
# Output:         List with mz, intensity, and metadata (filename, baseline).
# Usage:          Structures raw spectrum data for processing.

# read_spectrum_file(file_path)
# Purpose:        Reads spectrum data from .dat, .txt, or .csv files.
# Input:          file_path – Path to the file.
# Output:         Data frame with columns mz_value and intensity.
# Details:        Supports .dat, .txt (space-separated), and .csv formats; errors on unsupported formats.
# Usage:          Loads raw spectra in the Preprocessing panel.

# transform_intensity(sp, method='sqrt')
# Purpose:        Applies intensity transformation to a spectrum.
# Inputs:         sp – Spectrum object.
#                 method – Transformation method ('sqrt', 'log', or 'none').
# Output:         Spectrum object with transformed intensities.
# Details:        'sqrt' – Applies square root (sqrt(pmax(intensity, 0))).
#                 'log' – Applies log (log1p(pmax(intensity, 0))).
#                 Ensures non-negative intensities.
# Usage:          Optional preprocessing step in the Preprocessing panel.

# smooth_sgolay(y, polyOrder=3, windowSize=11)
# Purpose:        Applies Savitzky-Golay smoothing to intensity data.
# Inputs:         y – Intensity vector.
#                 polyOrder – Polynomial degree (default: 3).
#                 windowSize – Window size (default: 11, must be > polyOrder).
# Output:         Smoothed intensity vector.
# Details:        Uses signal::sgolayfilt for polynomial-based smoothing.
# Usage:          Optional smoothing method in the Preprocessing panel.

# baseline_snip(y, iterations=100)
# Purpose:        Estimates baseline using the SNIP algorithm.
# Inputs:         y – Intensity vector.
#                 iterations – Number of iterations (default: 100).
# Output:         Estimated baseline vector.
# Details:        Iteratively updates baseline by taking the minimum of current intensity and neighbor averages.
# Usage:          Baseline correction in the Preprocessing panel.

# remove_baseline(spectrum, baseline_vector)
# Purpose:        Subtracts baseline from spectrum intensities.
# Inputs:         spectrum – Spectrum object.
#                 baseline_vector – Estimated baseline vector.
# Output:         Spectrum object with baseline-corrected intensities.
# Details:        Ensures non-negative intensities using pmax.
# Usage:          Applied after baseline estimation in the Preprocessing panel.

# normalize_tic(spectrum)
# Purpose:        Normalizes intensities by Total Ion Current (TIC).
# Input:          Spectrum object.
# Output:         Normalized spectrum object.
# Details:        Divides intensities by their sum if non-zero.
# Usage:          Optional normalization in the Preprocessing panel.

# normalize_rmse(spectrum)
# Purpose:        Normalizes intensities by Root Mean Square Error.
# Input:          Spectrum object.
# Output:         Normalized spectrum object.
# Details:        Divides intensities by the square root of mean squared intensity.
# Usage:          Optional normalization in the Preprocessing panel.

# my_detect_peaks(spectrum, SNR=2, halfWindowSize=20)
# Purpose:        Detects peaks using the MAD method.
# Inputs:         spectrum – Spectrum object.
#                 SNR – Signal-to-noise ratio (default: 2).
#                 halfWindowSize – Half window size for peak detection (default: 20).
# Output:         Data frame with mz and intensity of detected peaks.
# Details:        Identifies local maxima.
#                 Filters peaks above SNR * mad(intensity).
# Usage:          Peak detection in the PeaksCatcher panel.

# find_bounds_valley(x, y, peak_x, min_distance=0.8, window=NULL)
# Purpose:        Finds peak boundaries using a valley-based method.
# Inputs:         x – m/z vector.
#                 y – Intensity vector.
#                 peak_x – m/z of peak apex.
#                 min_distance – Minimum distance for boundaries (default: 0.8).
#                 window – Optional m/z range (default: NULL).
# Output:         Indices vector of start and end boundaries.
# Details:        Uses derivative-based minima detection to locate valleys.
#                 Adjusts min_distance based on estimated FWHM.
#                 Ensures symmetric boundary detection within a specified window.
# Usage:          Determines peak integration ranges in the PeaksCatcher panel.

# calc_auc(data, start_x, end_x, baseline=0)
# Purpose:        Calculates the area under the curve (AUC) for a peak.
# Inputs:         data – Data frame with x (m/z) and y (intensity).
#                 start_x, end_x – Integration boundaries.
#                 baseline – Baseline intensity (default: 0).
# Output:         AUC value.
# Details:        Uses trapezoidal rule for integration.
# Usage:          Computes peak area in the PeaksCatcher panel.

# calc_peak_params(data, start_x, end_x, baseline=0)
# Purpose:        Calculates peak parameters (e.g., centroid, FWHM).
# Inputs:         data – Data frame with x (m/z) and y (intensity).
#                 start_x, end_x – Integration range.
#                 baseline – Baseline intensity (default: 0).
# Output:         List with peak_x, peak_y, centroid, FWHM, FWHM_left, FWHM_right.
# Details:        Computes peak apex, intensity-weighted centroid, and full width at half maximum.
# Usage:          Characterizes detected peaks in the PeaksCatcher panel.

# ------------------------------------------------------------------------------
