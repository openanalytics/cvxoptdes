AC_INIT([cvxoptdes],[1.0.1])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_PROG_CC

AC_PATH_PROG([CMAKE], [cmake], [],
  [$PATH:/Applications/CMake.app/Contents/bin])
if test -z "$CMAKE"; then
    AC_MSG_ERROR([no 'cmake' command found; install CMake or set CMAKE])
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
