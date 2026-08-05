AC_INIT([cvxoptdes],[1.0.0.9000])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_PROG_CC

if ! test -z "$CMAKE"; then
    echo "CMAKE set to: $CMAKE"
else 
    # Check if cmake is in the PATH
    CMAKE_PATH="$(command -v cmake)"
    if ! test -z "$CMAKE_PATH"; then
        CMAKE="$CMAKE_PATH"
        echo "CMAKE set to :$CMAKE"
    # If not, check the typical macOS CMake.app path
    elif [ -x "/Applications/CMake.app/Contents/bin/cmake" ]; then
        CMAKE="/Applications/CMake.app/Contents/bin/cmake"
        echo "CMAKE set to: $CMAKE"
    # If still not found, print error message
    else
        echo "no 'cmake' command found"
        exit 1
    fi
fi

if ! test -z "$MKLROOT"; then
    echo "MKLROOT set to: $MKLROOT"
    MKLROOT=${MKLROOT}
    AC_SUBST(MKLROOT)
fi

if ! test -z "$HIGHS_HOME"; then
    echo "HIGHS_HOME set to: $HIGHS_HOME"
    HIGHS_HOME=${HIGHS_HOME}
    AC_SUBST(HIGHS_HOME)
else 
    echo "HIGHS_HOME undefined, not using external HiGHS installation"
fi

AC_SUBST(CMAKE)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
