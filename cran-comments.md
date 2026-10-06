## CRAN package version 1.0.0

* System requirements: GNU make, CMake

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