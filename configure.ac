# Process this file with autoconf to produce a configure script. 
AC_INIT(BHC, 1.2.0) 
AC_CONFIG_SRCDIR([src/bhc.cpp])
AC_PREREQ([2.62])

# find R home and set CC/CFLAGS 
: ${R_HOME=`R RHOME`} 
if test -z "${R_HOME}"; then 
  echo "could not determine R_HOME" 
  exit 1 
fi 
RBIN="${R_HOME}/bin/R" 
CC=`"${RBIN}" CMD config CC`
CFLAGS=`"${RBIN}" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LIBS="${PKG_LIBS}" 

# Check for programs
AC_PROG_CC 
# Check for OpenMP 
AC_OPENMP 

# since some systems have broken OMP libraries 
# we also check that the actual package will work 
ac_pkg_openmp=no 
if test -n "${OPENMP_CFLAGS}"; then 
  AC_MSG_CHECKING([whether OpenMP will work in this package]) 
  AC_LANG_CONFTEST( 
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads(); ]])]) 
  PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes 
  AC_MSG_RESULT([${ac_pkg_openmp}]) 
fi 

# if ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no" 
if test "${ac_pkg_openmp}" = no; then 
  OPENMP_CFLAGS=''
  AC_MSG_WARN([No OpenMP support. If using GCC, upgrade to >= 4.2])
fi 

AC_SUBST(OPENMP_CFLAGS) 

AC_CONFIG_FILES([src/Makevars]) 
AC_OUTPUT  
