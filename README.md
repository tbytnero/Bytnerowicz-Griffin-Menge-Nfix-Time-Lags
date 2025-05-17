# Bytnerowicz-Griffin-Menge-Nfix-Time-Lags

R scripts and data to recreate analyses and figures for:

Bytnerowicz et al. "Time lags in the regulation of symbiotic nitrogen fixation"

There is one R script.

The R script recreates the analyses and figures for the main text and for the supplement.

There are four .csv files of data:

The first is nitrogen fixation data (normalized by whole-symbiosis respiration): snf_time_lags.csv

The second is the photosynthesis and respiration data: co2_time_lags.csv

The third is the rinse test data: rinse_test.csv

The fourth is the raw nitrogen fixation data (not normalized by whole-symbiosis respiration): raw_Nfix.csv

The nitrogen fixation rates in the uploaded data (first and fourth file) still need to be converted to units of N2 fixed from acetylene reduced. This is accomplished in the R script using the values in Table S3.
