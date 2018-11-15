AC_INIT([wgms3d], [2.0])
AC_CONFIG_SRCDIR([wgms3d.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER(config.h)

##############################################################################
# Check for compilers (and MPI, if requested)
#

AC_ARG_WITH([mpi], AS_HELP_STRING([--with-mpi], [Build with MPI (default: disabled)]))
AS_IF([test "x$with_mpi" = "xyes"], [F77="mpif77" CXX="mpic++"])
AC_PROG_F77
AC_PROG_CXX
AS_IF([test "x$with_mpi" = "xyes"], [ LDFLAGS="$LDFLAGS -lmpi" CPPFLAGS="$CPPFLAGS -I/usr/include/mpi"
                                      AC_DEFINE(HAVE_MPI,1,[Defined if we're compiling against MPI.]) ])
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS

AC_LANG(C++)
sinclude(ax_cxx_compile_stdcxx_11.m4)
AX_CXX_COMPILE_STDCXX_11

AC_CHECK_HEADER([boost/throw_exception.hpp], [], [AC_MSG_ERROR(Need Boost C++ library headers.)])

##############################################################################
# Check for BLAS and LAPACK
#

AC_LANG(C)

sinclude(acx_blas.m4)
sinclude(acx_lapack.m4)

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##############################################################################
# Check to see if calling Fortran functions (in particular, the BLAS
# complex dot product functions) works from C:
# (adapted from mpb-1.4.2)

AC_MSG_CHECKING([whether calling BLAS zdotc works])
AC_F77_FUNC(zdotc)
AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
typedef struct { double re, im; } cmplx;
#ifdef __cplusplus
extern "C"
#endif
extern cmplx $zdotc(int *, cmplx *, int *, cmplx *, int *);
int main(void) {
     cmplx x[2] = { {1,0}, {0,1} }, y[2] = { {1,0}, {1,0} };
     int n = 2, inc1 = 1, inc2 = 1;
     cmplx result = { 0, 0 };
     result = $zdotc(&n, x, &inc1, y, &inc2);
     if (result.re != 1 || result.im != -1) return EXIT_FAILURE;
     else return EXIT_SUCCESS;
}
], ok="yes", ok="no", ok="no")
AC_MSG_RESULT($ok)
if test "$ok" = "no"; then
	AC_MSG_ERROR(Your BLAS doesn't work.)
fi

##############################################################################
# Check for SuperLU
#

AC_ARG_WITH(
 superlu,
 AS_HELP_STRING([--with-superlu=DIR], [Build with SuperLU library (default: disabled)]),
 [LDFLAGS="$LDFLAGS -L${withval}/lib" CPPFLAGS="$CPPFLAGS -I${withval}/include"],
 [with_superlu=no]
)
AS_IF([test "x$with_superlu" != "xno"], [
 AC_CHECK_LIB([superlu], [dgstrf], [have_superlu=yes; LIBS="-lsuperlu $LIBS"; AC_DEFINE(HAVE_SUPERLU,1,[Defined if we can link to SUPERLU.])])
 AS_IF([test "x$have_superlu" != "xyes"], [ AC_MSG_ERROR(Please try --with-superlu={prefix}.)])
])

##############################################################################
# Check for ARPACK
#

AC_ARG_WITH(
 arpack,
 AS_HELP_STRING([--with-arpack=DIR], [Build with ARPACK library (default: disabled)]),
 [LDFLAGS="$LDFLAGS -L${withval}/lib" ],
 [with_arpack=no]
)
AS_IF([test "x$with_arpack" != "xno"], [
 # First make sure that SuperLU is working:
 AS_IF([test "x$have_superlu" != "xyes"], [ AC_MSG_ERROR(Can't use ARPACK without SuperLU.) ])
 AC_F77_FUNC(dnaupd)
 AC_CHECK_LIB([arpack], [$dnaupd], [have_arpack=yes; LIBS="-larpack $LIBS"; AC_DEFINE(HAVE_ARPACK,1,[Defined if we can link to ARPACK.])])
 AS_IF([test "x$have_arpack" != "xyes"], [ AC_MSG_ERROR(Please try --with-arpack={prefix}.)])
])

AM_CONDITIONAL([HAS_ARPACK], [test "x$have_arpack" = "xyes"])

##############################################################################
# Check for PETSc
#

AC_LANG(C++)

AC_ARG_WITH(petsc, AS_HELP_STRING([--with-petsc], [Build with PETSc library (default: disabled)]))
AS_IF([test "x$with_petsc" = "xyes"], [
 AS_IF([test "x$with_mpi" != "xyes"], [ AC_MSG_ERROR(Can't use PETSc without MPI.) ])
 AS_IF([test "x$with_petsc" != "xno"], [
 	     petsclib=`(cd $PETSC_DIR; make getlinklibs)`;
 	     petscinclude=`(cd $PETSC_DIR; make getincludedirs)`;
	     LDFLAGS="$LDFLAGS $petsclib"
	     CPPFLAGS="$CPPFLAGS $petscinclude"])
 AC_CHECK_HEADER([petsc.h],
  [AC_CHECK_LIB([petsc], [MatCreate], [have_petsc=yes; LIBS="-lpetsc $LIBS"; AC_DEFINE(HAVE_PETSC,1,[Defined if we can link to PETSc.])])])
 AS_IF([test "x$have_petsc" != "xyes"], [ AC_MSG_ERROR(Please set PETSC_DIR and PETSC_ARCH environment variables.)])
])

##############################################################################
# Check for SLEPc
#

AC_ARG_WITH(slepc, AS_HELP_STRING([--with-slepc], [Build with SLEPc library (default: disabled)]))
AS_IF([test "x$with_slepc" = "xyes"], [
 AS_IF([test "x$have_petsc" != "xyes"], [AC_MSG_ERROR(Can't use SLEPc without PETSc.)])
 AS_IF([test "x$with_slepc" != "xno"], [
 	     #slepclib=`(cd $SLEPC_DIR; make getlinklibs)`;
 	     #slepcinclude=`(cd $SLEPC_DIR; make getincludedirs)`;
 	     LDFLAGS="$LDFLAGS -L${SLEPC_DIR}/${PETSC_ARCH}/lib -lpetsc"
	     CPPFLAGS="$CPPFLAGS -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include"])
 AC_CHECK_HEADER([slepceps.h],
  [AC_CHECK_LIB([slepc], [EPSCreate], [have_slepc=yes; LIBS="-lslepc $LIBS"; AC_DEFINE(HAVE_SLEPC,1,[Defined to 1 if we can link to SLEPc.])])])
 AS_IF([test "x$have_slepc" != "xyes"], [ AC_MSG_ERROR(Please set SLEPC_DIR environment variable.)])
])
AM_CONDITIONAL([HAS_SLEPC], [test "x$with_slepc" = "xyes"])

##############################################################################
# Check for collisions
#

AS_IF([test "x$have_superlu" = "xyes" -a "x$have_arpack" != "xyes"],
      [AC_MSG_ERROR(Having SuperLU without ARPACK doesn't make much sense.)])

AS_IF([test "x$have_petsc" = "xyes" -a "x$have_slepc" != "xyes"],
      [AC_MSG_ERROR(Having PETSc without SLEPc doesn't make much sense.)])

AS_IF([test "x$with_slepc" = "xyes" -a "x$have_arpack" = "xyes"],
      [AC_MSG_ERROR(Can't compile SuperLU/ARPACK and SLEPc solvers into wgms3d simultaneously.)])

AS_IF([test "x$have_slepc" != "xyes" -a "x$have_arpack" != "xyes"],
      [AC_MSG_ERROR(No solver selected.)])

##############################################################################
# Build Makefile
#

AC_CONFIG_FILES([Makefile
		 doc/Makefile
		 matlab/Makefile
		 tests/Makefile
		 tests/all_field_components/Makefile
		 tests/complex_modes/Makefile
		 tests/fiber_convergence/Makefile
		 tests/lossy_materials/Makefile
		 tests/semivectorial/Makefile
		 tests/disk_resonator/Makefile
		 tests/silicon_strip_waveguide/Makefile
		 tests/tm2te_leakage/Makefile])
AC_OUTPUT
