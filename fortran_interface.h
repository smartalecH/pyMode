
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _FORTRAN_INTERFACE_H
#define _FORTRAN_INTERFACE_H

#include "config.h"

#include "complex_functions.h"

extern "C"
{
    /* BLAS/LAPACK functions */
    extern void F77_FUNC(zgels,ZGELS) (char *TRANS, int *M, int *N, int *NRHS, std::complex<double> *A, int *LDA,
				       std::complex<double> *B, int *LDB, std::complex<double> *WORK, int *LWORK, int *INFO);
    extern void F77_FUNC(dcopy,DCOPY) (int *n, double *x, int *incx, double *y, int *incy);
    extern void F77_FUNC(zcopy,ZCOPY) (int *n, std::complex<double> *x, int *incx, std::complex<double> *y, int *incy);
    extern void F77_FUNC(dscal,DSCAL) (int *n, double *alpha, double *x, int *incx);
    extern void F77_FUNC(zscal,ZSCAL) (int *n, std::complex<double> *alpha, std::complex<double> *x, int *incx);
    extern void F77_FUNC(daxpy,DAXPY) (int *n, double *alpha, double *x,
				       int *incx, double *y, int *incy);
    extern void F77_FUNC(zaxpy,ZAXPY) (int *n, std::complex<double> *alpha, std::complex<double> *x,
				       int *incx, std::complex<double> *y, int *incy);
    extern void F77_FUNC(dgemm,DGEMM) (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
				       double *ALPHA, double *A, int *LDA, double *B, int *LDB,
				       double *BETA, double *C, int *LDC);
    extern void F77_FUNC(zgemm,ZGEMM) (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
				       std::complex<double> *ALPHA, std::complex<double> *A, int *LDA,
				       std::complex<double> *B, int *LDB,
				       std::complex<double> *BETA, std::complex<double> *C, int *LDC);
}

inline void COPY (int n, double *x, int incx, double *y, int incy) {
    F77_FUNC(dcopy,DCOPY) (&n, x, &incx, y, &incy); }
inline void COPY (int n, std::complex<double> *x, int incx, std::complex<double> *y, int incy) {
    F77_FUNC(zcopy,ZCOPY) (&n, x, &incx, y, &incy); }
inline void SCAL (int n, double alpha, double *x, int incx) {
    F77_FUNC(dscal,DSCAL) (&n, &alpha, x, &incx); }
inline void SCAL (int n, std::complex<double> alpha, std::complex<double> *x, int incx) {
    F77_FUNC(zscal,ZSCAL) (&n, &alpha, x, &incx); }
inline void GELS (char *TRANS, int *M, int *N, int *NRHS, std::complex<double> *A, int *LDA,
		  std::complex<double> *B, int *LDB, std::complex<double> *WORK, int *LWORK, int *INFO) {
    F77_FUNC(zgels,ZGELS) (TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO); }
inline void AXPY (int n, double alpha, double *x,
		  int incx, double *y, int incy) {
    F77_FUNC(daxpy,DAXPY) (&n, &alpha, x, &incx, y, &incy); }
inline void AXPY (int n, std::complex<double> alpha, std::complex<double> *x,
		  int incx, std::complex<double> *y, int incy) {
    F77_FUNC(zaxpy,ZAXPY) (&n, &alpha, x, &incx, y, &incy); }
inline void GEMM (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
		  double *ALPHA, double *A, int *LDA, double *B, int *LDB,
		  double *BETA, double *C, int *LDC) {
    F77_FUNC(dgemm,DGEMM) (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC); }
inline void GEMM (char *TRANSA, char *TRANSB, int *M, int *N, int *K,
		  std::complex<double> *ALPHA, std::complex<double> *A, int *LDA,
		  std::complex<double> *B, int *LDB,
		  std::complex<double> *BETA, std::complex<double> *C, int *LDC) {
    F77_FUNC(zgemm,ZGEMM) (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC); }

template <class T>
void dump_fortran_matrix (T *A,
			  int lda,
			  int nrows,
			  int ncols)
{
    int i, j;

    for(i = 0; i < nrows; i++) {
	for(j = 0; j < ncols; j++) {
	    //	    printf("% 2.1e+i% 2.1e ", A[i+j*lda].real(), A[i+j*lda].imag());
	    if(MABS(A[i+j*lda]) < 1*1e-12) {
		printf("          ");
	    } else {
		printf("% 2.2e ", A[i+j*lda]);
	    }
	}
	printf("\n");
    }
}

template <class T>
void dump_fortran_matrix (std::complex<T> *A,
			  int lda,
			  int nrows,
			  int ncols)
{
    int i, j;

    for(i = 0; i < nrows; i++) {
	for(j = 0; j < ncols; j++) {
	    printf("% 2.1e+i% 2.1e ", A[i+j*lda].real(), A[i+j*lda].imag());
	}
	printf("\n");
    }
}

template <class T>
void dump_fortran_matrix (T *A,
			  int nrows,
			  int ncols)
{
    dump_fortran_matrix(A, nrows, nrows, ncols);
}

template <class T>
void matrix_copy (T *dst,
		  int ldst,
		  T *src,
		  int lsrc,
		  int nrows,
		  int ncols)
{
    int i, j;

    for(i = 0; i < nrows; i++) {
	for(j = 0; j < ncols; j++) {
	    dst[i + j*ldst] = src[i + j*lsrc];
	}
    }
}

#endif // _FORTRAN_INTERFACE_H

