## Test environments
* local OS X install, R 3.4.0
* ubuntu 12.04 (on travis-ci), R 3.4.0
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs

## Resubmission

This is a resubmission. In this version I have 
- Fixed a bug in output changepoints from simulateAR1
- Added another output to estimateSpikes
  * Now output cost F(s) at each data point s = 1, ..., T
  

