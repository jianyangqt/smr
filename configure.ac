#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([smr], [1.0], [])
AC_PREREQ([2.69])
AC_CONFIG_SRCDIR([CommFunc.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lomp':
AC_CHECK_LIB([omp], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset mkdir pow setenv sqrt strerror])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
