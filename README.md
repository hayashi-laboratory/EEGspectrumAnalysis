# EEGspectrumAnalysis

[![DOI](https://zenodo.org/badge/334886963.svg)](https://zenodo.org/badge/latestdoi/334886963)

EEG power spectrum analysis from SleepSign
FFTspec5: calculate the pure Wake, NREM ,REM (REM>40secs)
FFTspec6: calculate all pure Wake, NREM , REM (all length)

# Use FFTspec5Run to execute FFTspec5

# Use FFTspec6Run to execute FFTspec6

FFTspec5/6: 
(1) Run text (.txt) files, which will be very fast.

​	Method of exporting FFT .txt file from SleepSign: 

​	Analysis → FFT parameter → EpochFF/AverageFFT (Average:4, Display range: 0.0 Hz-30 Hz) → ok

​	Analysis → FFT Text Output → Continuous FFT → epoch XXX (start)-epoch XXX (end) 

​	export: File, Path (export to the same folder) → export as .txt file (a lot faster to read by MATLAB) → ✔︎comma(csv) → start

(2) Can use in AE(artifact exclusion) file, since M will be viewed as NaN
