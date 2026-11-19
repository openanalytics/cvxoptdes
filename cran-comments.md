## CRAN package version 1.0.1

* System requirements: GNU make, CMake
* Fixes CRAN check issues:
  - LTO/gcc-asan: fixed SCS compilation warnings
  - macos-x86_64/macos-arm64: fix issue locating CMake (macos) in configure

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
