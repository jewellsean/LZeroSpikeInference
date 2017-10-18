devtools::build_win()## Test environments
* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs

## Resubmission

This is a resubmission. In this version I have 
- Fixed a bug in output changepoints from simulateAR1
- Edited simulateAR1 so that c_1 = 0
- Added another output to estimateSpikes
  * Now output cost F(s) at each data point s = 1, ..., T
- Updated the optimization problem to include hardThresholding for estimated calcium concentration
- Fixed a bug in the cross validation procedure
  

