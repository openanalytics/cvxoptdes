## CRAN package version 1.0.0

* System requirements: GNU make, CMake
* Fixes initial submission reducing tarball size to <10MB

## Test environments

* ubuntu gcc/clang R-release, R-next, R-devel (local install, rhub)
* fedora gcc R-devel (rhub)
* macos gcc/clang R-devel (rhub)
* windows gcc R-oldrel, R-release, R-devel (r-winbuilder)

## Compiled code checks

* ubuntu-rchk R-devel (docker)
* ubuntu gcc R-release --use-valgrind (local install) 
* macos clang ASAN/UBSAN R-devel (rhub)
* ubuntu clang ASAN/UBSAN R-devel (rhub)

## R CMD check results

1 NOTE (ubuntu gcc/clang): 
* Running R code in ‘test_cvxoptdes.R’ had CPU time 2.8 times elapsed time
