# EEGPowerCalc
 EEG power analysis from SleepSign

FFTspec6: 



(1) Run text (.txt) files, which will be very fast.



​	Method of exporting FFT .txt file from SleepSign: 

​	Analysis → FFT parameter → EpochFF/AverageFFT (Average:4, Display range: 0.0 Hz-30 Hz) → ok

​	Analysis → FFT Text Output → Continuous FFT → epoch XXX (start)-epoch XXX (end) 

​	export: File, Path (export to the same folder) → export as .txt file (a lot faster to read by MATLAB) → ✔︎comma(csv) → start



(2) Can use in AE(artifact exclusion) file, since M will be viewed as NaN