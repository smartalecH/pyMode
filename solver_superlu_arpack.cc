
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2014  Michael Krause <m.krause@tu-harburg.de>

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

#include "config.h"

#include <algorithm>
#include <complex>

#include <slu_ddefs.h>
#include <slu_util.h>

#include "fortran_interface.h"
#include "solver.h"
#include "complex_functions.h"
#include "make_unique.h"

extern "C" {

    /* The following works on Linux, and on HPUX in 32-bit and 64-bit modes */
    typedef int logical;

    /* ARPACK functions */
    extern void
    F77_FUNC(znaupd,ZNAUPD) (int *ido, char *bmat, int *n, char *which, int *nev,
			     double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v,
			     int *ldv, int *iparam, int *ipntr, std::complex<double> *workd,
			     std::complex<double> *workl, int *lworkl, double *rwork, int *info);
    extern void
    F77_FUNC(dnaupd,DNAUPD) (int *ido, char *bmat, int *n, char *which, int *nev,
			     double *tol, double *resid, int *ncv, double *v,
			     int *ldv, int *iparam, int *ipntr, double *workd,
			     double *workl, int *lworkl, int *info);
    extern void
    F77_FUNC(zneupd,ZNEUPD) (logical *rvec, char *howmny, logical *select, std::complex<double> *d,
			     std::complex<double> *z, int *ldz, std::complex<double> *sigma, std::complex<double> *workev,
			     char *bmat, int *n, char *which, int *nev, double *tol,
			     std::complex<double> *resid, int *ncv, std::complex<double> *v, int *ldv,
			     int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl,
			     int *lworkl, double *rwork, int *info);
    extern void
    F77_FUNC(dneupd,DNEUPD) (logical *rvec, char *howmny, logical *select, double *dr, double *di,
			     double *z, int *ldz, double *sigmar, double *sigmai, double *workev,
			     char *bmat, int *n, char *which, int *nev, double *tol, double *resid,
			     int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd,
			     double *workl, int *lworkl, int *info);

    /* SuperLU functions */
    /* We want to use both the real and complex SuperLU routines, but we
       can't include slu_ddefs.h and slu_zdefs.h at the same time, since
       they contain colliding definitions.

       Therefore we copy here the definitions from slu_zdefs.h...

       This is for SuperLU Version 4.0. */

    extern void    zgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
			   SuperMatrix *, SuperLUStat_t*, int *);
    
    extern void    zgstrf (superlu_options_t*, SuperMatrix*,
			   int, int, int*, void *, int, int *, int *, 
			   SuperMatrix *, SuperMatrix *, SuperLUStat_t*, int *);

    extern void
    zCreate_CompCol_Matrix(SuperMatrix *, int, int, int, std::complex<double> *,
			   int *, int *, Stype_t, Dtype_t, Mtype_t);

    extern int     
    zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
}

namespace
{

    void GSTRS (trans_t a, SuperMatrix *b, SuperMatrix *c, int *d, int *e,
		SuperMatrix *B, SuperLUStat_t *f, int *g)
    {
	if(B->Dtype == SLU_D) {
	    dgstrs(a,b,c,d,e,B,f,g);
	} else if(B->Dtype == SLU_Z) {
	    zgstrs(a,b,c,d,e,B,f,g);
	}
    }

    void GSTRF (superlu_options_t *a, SuperMatrix *AC,
		int d, int e, int *f, void *g, int h, int *i, int *j, 
		SuperMatrix *k, SuperMatrix *l, SuperLUStat_t *m, int *n)
    {
	if(AC->Dtype == SLU_D) {
	    dgstrf(a,AC,d,e,f,g,h,i,j,k,l,m,n);
	} else if(AC->Dtype == SLU_Z) {
	    zgstrf(a,AC,d,e,f,g,h,i,j,k,l,m,n);
	}
    }

    int QuerySpace (SuperMatrix *L, SuperMatrix *U, mem_usage_t *m)
    {
	if(L->Dtype == SLU_D) {
	    dQuerySpace(L, U, m);
	} else if(L->Dtype == SLU_Z) {
	    zQuerySpace(L, U, m);
	}
	return 0;
    }

    inline void
    NAUPD (int *ido, char *bmat, int *n, char *which, int *nev,
	   double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v,
	   int *ldv, int *iparam, int *ipntr, std::complex<double> *workd,
	   std::complex<double> *workl, int *lworkl, double *rwork, int *info)
    {
	F77_FUNC(znaupd,ZNAUPD) (ido, bmat, n, which, nev, tol, resid, ncv, v,
				 ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
    }

    inline void
    NAUPD (int *ido, char *bmat, int *n, char *which, int *nev,
	   double *tol, double *resid, int *ncv, double *v,
	   int *ldv, int *iparam, int *ipntr, double *workd,
	   double *workl, int *lworkl, double *rwork, int *info)
    {
	F77_FUNC(dnaupd,DNAUPD) (ido, bmat, n, which, nev, tol, resid, ncv, v,
				 ldv, iparam, ipntr, workd, workl, lworkl, info);
    }

    inline void
    NEUPD (logical *rvec, char *howmny, logical *select, std::complex<double> *d,
	   std::complex<double> *z, int *ldz, std::complex<double> sigma, std::complex<double> *workev,
	   char *bmat, int *n, char *which, int *nev, double *tol,
	   std::complex<double> *resid, int *ncv, std::complex<double> *v, int *ldv,
	   int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl,
	   int *lworkl, double *rwork, int *info)
    {
	F77_FUNC(zneupd,ZNEUPD) (rvec, howmny, select, d, z, ldz, &sigma,
				 workev, bmat, n, which, nev, tol, resid, ncv,
				 v, ldv, iparam, ipntr, workd, workl, lworkl,
				 rwork, info);
    }

    inline void
    NEUPD (logical *rvec, char *howmny, logical *select, std::complex<double> *d,
	   double *z, int *ldz, double sigma, double *workev,
	   char *bmat, int *n, char *which, int *nev, double *tol,
	   double *resid, int *ncv, double *v, int *ldv,
	   int *iparam, int *ipntr, double *workd, double *workl,
	   int *lworkl, double *rwork, int *info)
    {
	int i;
	int nevcopy = *nev;
    
	double *dr = new double[nevcopy+1];
	double *di = new double[nevcopy+1];
	    
	double sigmar = sigma;
	double sigmai = 0.0;
	    
	F77_FUNC(dneupd,DNEUPD) (rvec, howmny, select, dr, di, z, ldz, &sigmar, &sigmai,
				 workev, bmat, n, which, nev, tol, resid, ncv,
				 v, ldv, iparam, ipntr, workd, workl, lworkl, info);
	    
	for(i = 0; i < nevcopy+1; i++) {
	    d[i] = std::complex<double>(dr[i],di[i]);
	}
	    
	delete[] dr;
	delete[] di;
    }

    inline void
    Create_CompCol_Matrix (SuperMatrix *a, int b, int c, int d, double *e,
			   int *f, int *g, Stype_t h, Mtype_t j)
    {
	dCreate_CompCol_Matrix (a, b, c, d, e, f, g, h, SLU_D, j);
    }
    
    inline void
    Create_CompCol_Matrix (SuperMatrix *a, int b, int c, int d, std::complex<double> *e,
			   int *f, int *g, Stype_t h, Mtype_t j)
    {
	zCreate_CompCol_Matrix (a, b, c, d, reinterpret_cast <std::complex<double>*> (e), f, g, h, SLU_Z, j);
    }

    inline Dtype_t
    get_slu_type (double *)
    {
	return SLU_D;
    }
    inline Dtype_t
    get_slu_type (std::complex<double> *)
    {
	return SLU_Z;
    }

    // RAII wrapper
    class LUDecomposition
    {
    public:
	SuperMatrix L;
	SuperMatrix U;
	int *perm_c;
	int *perm_r;

	LUDecomposition()
	    : perm_c(nullptr),
	      perm_r(nullptr)
	{
	}

	~LUDecomposition() {
	    Destroy_SuperNode_Matrix(&L);
	    Destroy_CompCol_Matrix(&U);
	    delete[] perm_c;
	    delete[] perm_r;
	}
    };

    /* superlu(): convert the sparse FD system matrix (always in
     * complex double format) A to SuperLU format (in the format
     * specified by T). */
    template <class T>
    bool superlu (sparse_matrix<std::complex<double> > &A,
		  LUDecomposition &result,
		  T /* tag */)
    {
	SuperLUStat_t stat;
	superlu_options_t options;
	SuperMatrix A_slu;
	int permc_spec;
	int nnz = A.length;
	T *a = NULL;
	int *asub = NULL, *xa = NULL;
	int ap = 0;
	int *etree = NULL;
	SuperMatrix AC;
	int panel_size;
	int relax;
	int lwork = 0;
	int info;
	unsigned int i, j;
	mem_usage_t mem_usage;

	a = new T[nnz];
	asub = new int[nnz];
	xa = new int[A.n+1];
	etree = new int[A.n];

	panel_size = sp_ienv(1);
	relax = sp_ienv(2);

	/* Convert system matrix to SuperLU format */
	A.order(2);
	ap = 0;
	for(j = 0; j < A.n; j++) {
	    xa[j] = ap;
	    for(i = A.indextable[j]; i < A.indextable[j+1]; i++) {
		a[ap] = wgms3d::maybeConvertComplexToReal(A.entries[i].v, T());
		asub[ap] = A.entries[i].i;
		ap++;
	    }
	}
	xa[j] = ap;

	set_default_options(&options);
	options.ColPerm = COLAMD;

	Create_CompCol_Matrix(&A_slu, A.n, A.n, ap, a, asub, xa, SLU_NC, SLU_GE);

	result.perm_c = new int[A.n];
	permc_spec = options.ColPerm;
	if (permc_spec != MY_PERMC && options.Fact == DOFACT)
	    get_perm_c(permc_spec, &A_slu, result.perm_c);

	StatInit(&stat);

	sp_preorder(&options, &A_slu, result.perm_c, etree, &AC);

	result.perm_r = new int[A.n];

	GSTRF(&options, &AC, relax, panel_size,
	      etree, NULL, lwork, result.perm_c, result.perm_r, &result.L, &result.U, &stat, &info);
	if(info) {
	    std::cerr << "SuperLU Xgstrf failed with INFO = " << info << std::endl;
	    exit(1);
	}

	Destroy_SuperMatrix_Store(&A_slu);
	Destroy_CompCol_Permuted(&AC);

	QuerySpace(&result.L, &result.U, &mem_usage);
	std::cout << " (~" << (int)(mem_usage.total_needed/(1<<20)) << "MB)" << std::endl;

	delete[] etree;
	delete[] xa;
	delete[] asub;
	delete[] a;

	return true;
    }

    template <class T>
    bool eigensolve (LUDecomposition &lu,
		     double sigma,
		     std::unique_ptr< T[] > &pEvec,                    /* Store pointer to eigenvector array here (on success) */
		     std::unique_ptr< std::complex<double>[] > &pEval, /* Store pointer to eigenvalue array here (on success) */
		     int nev) // in eigs.m: 'k'
    {
	int ido = 0;
	char bmat = 'I';
	int n = lu.L.nrow; // dimension of system matrix
	char which[] = "LM";
	double tol = 0.0;
	T *resid = new T[2*n]; // n
	int ncv = std::min(std::max(2*nev+1,20),n); // from Matlab's eigs.m, there: 'p'
	int ldv = n;
	int iparam[11];
	int ipntr[15];
	T *workd = NULL; 
	int lworkl = 3*ncv*ncv + 6*ncv;
	T *workl = new T[lworkl];
	int info = 0;
	bool rc = false;
	T sigma2 = sigma;

	logical rvec = 1;
	char howmny = 'A';
	logical *select = new logical[ncv];
	T *workev = new T[3*ncv];
	double *rwork = new double[ncv];

	/* SuperLU variables */
	SuperLUStat_t stat;
	int slu_info;

	if(ncv < nev + 2 || ncv > n) {
	    std::cerr << "Too many eigenvalues requested." << std::endl;
	    exit(1);
	}

	memset(iparam, 0, sizeof(iparam));
	iparam[0] = 1;
	iparam[2] = 10000;
	iparam[6] = 3; // shift-invert mode, M = id

	std::cout << "Eigensolving using ARPACK (nev=" << nev << ", ncv=" << ncv << ")..." << std::endl;

	workd = new T[3*n];
	pEvec.reset(new T[ldv*ncv]);
	pEval.reset(new std::complex<double>[nev+1]);

	StatInit(&stat);

	while(1) {
	    NAUPD(&ido,
		  &bmat, &n, which, &nev, &tol, resid, &ncv, pEvec.get(), &ldv,
		  iparam, ipntr, workd, workl, &lworkl, rwork, &info);

	    if(ido == -1 || ido == 1) {
		/* Solve linear system with (A - \sigma I) */
		T *rhs = workd + ipntr[0] - 1;
		T *dst = workd + ipntr[1] - 1;

		trans_t trans = NOTRANS;
		NCformat Bstore;
		SuperMatrix B;

		Bstore.nnz = n;
		Bstore.nzval = dst;
		Bstore.rowind = NULL;
		Bstore.colptr = NULL;
		B.Stype = SLU_DN;
		B.Dtype = get_slu_type(rhs);
		B.Mtype = SLU_GE;
		B.nrow = n;
		B.ncol = 1;
		B.Store = &Bstore;

		memcpy(dst, rhs, n * sizeof(T));
		GSTRS(trans, &lu.L, &lu.U, lu.perm_c, lu.perm_r, &B, &stat, &slu_info);
		if(slu_info) {
		    std::cerr << "SuperLU Xgstrs failed with INFO = "
			      << slu_info << std::endl;
		    exit(1);
		}
	    } else {
		break;
	    }
	}

	if(ido != 99 || info != 0) {
	    std::cerr << "ARPACK Xnaupd finished with IDO = "
		      << ido << ", INFO = " << info << std::endl;
	    exit(1);
	}

	std::cout << "Eigencalculation finished successfully (niter="
		  << iparam[2] << ", nconv=" << iparam[4] << std::endl;

	NEUPD(&rvec, &howmny, select, pEval.get(), pEvec.get(), &ldv,
	      sigma2, workev,
	      &bmat, &n, which, &nev, &tol, resid, &ncv, pEvec.get(), &ldv,
	      iparam, ipntr, workd, workl, &lworkl, rwork, &info);
	if(info != 0) {
	    std::cout << "Error with dneupd, INFO = " << info << std::endl;
	    exit(1);
	}

	rc = true;

	delete[] resid;
	delete[] workl;
	delete[] select;
	delete[] workev;
	delete[] rwork;
	delete[] workd;

	return rc;
    }

    /* Copy complex-conjugate eigenvector returned by DNEUPD to a
     * complex array, optionally conjugating the vector in the
     * process. */
    static void
    cpfield (std::complex<double> *to,
	     double *dneupd_vector,
	     int realpartoffset,
	     int imagpartoffset,
	     int imagpartmult,
	     int fcsize)
    {
	int i;

	if(!dneupd_vector)
	    return;

	for(i = 0; i < fcsize; i++) {
	    to[i] = std::complex<double>(dneupd_vector[i+realpartoffset],
				    imagpartmult*dneupd_vector[i+imagpartoffset]);
	}
    }

    static void
    get_complex_ht_sub (int modetype,
			std::complex<double> *Ht,
			double *evec,
			int evecsize,
			int n)
    {
	switch(modetype) {
	case 0:
	    /* first of complex-conjugate pair:
	     * Hxy0[0:fcsize-1] is real part,
	     * Hxy0[evecsize:evecsize+fcsize-1] is imaginary
	     * part */
	    cpfield(Ht, evec, 0, evecsize, 1, n);
	    break;
	case 1:
	    /* second of pair:
	     * Hxy0[-evecsize:-evecsize+fcsize-1] is real
	     * part, Hxy0[0:fcsize-1] is imaginary part (to be
	     * conjugated!) */
	    cpfield(Ht, evec, -evecsize, 0, -1, n);
	    break;
	case -1:
	    /* Real eigenvalue. Just copy over the field. There's
	     * a real part only. */
	    int j;
	    for(j = 0; j < n; j++) {
		Ht[j] = evec[j];
	    }
	    break;
	default:
	    assert(1 == 0);
	    break;
	}
    }

    class ArpackSolverResult : public ISolverResult
    {
    private:
	const int n;
	const size_t num_eigenpairs;
	const bool is_matrix_real;

	LUDecomposition lu;
	std::unique_ptr< std::complex<double>[] > pEval;
	std::unique_ptr< std::complex<double>[] > pEvecComplex;
	std::unique_ptr< double[]               > pEvecReal;
	std::vector<int> modetypes;

    public:
	ArpackSolverResult (
	    std::unique_ptr< sparse_matrix<std::complex<double>> > matrix,
	    bool is_matrix_real,
	    double sigma,
	    int num_eigenvalues)
	    : n(matrix->n),
	      num_eigenpairs(num_eigenvalues),
	      is_matrix_real(is_matrix_real)
	{
	    matrix->add_constdiag(-sigma);

	    bool result;

	    /* Use SuperLU/ARPACK functions optimized for real matrices if
	     * possible (faster, less memory consumption). */
    
	    if(is_matrix_real) {
		std::cout << "Factorizing FD matrix using (real) SuperLU..." << std::flush;
		result = superlu(*matrix, lu, double());
	    } else {
		std::cout << "Factorizing FD matrix using (complex) SuperLU..." << std::flush;
		result = superlu(*matrix, lu, std::complex<double>());
	    }
	    matrix.reset();
	    if(!result) {
		BOOST_THROW_EXCEPTION(std::runtime_error("Couldn't LU-factorize system matrix."));
	    }

	    if(is_matrix_real) {
		result = eigensolve(lu, sigma, pEvecReal,    pEval, num_eigenvalues);
	    } else {
		result = eigensolve(lu, sigma, pEvecComplex, pEval, num_eigenvalues);
	    }
	    if(!result) {
		BOOST_THROW_EXCEPTION(std::runtime_error("Couldn't perform eigensolve."));
	    }

	    int modetype = 0;
	    for(int i = 0; i < num_eigenvalues; i++) {
		std::complex<double> beta2 = pEval[i]; /* This is beta^2 */

		/* the following is meaningful only for the real calculation
		 * modes. this is to categorize the modes as:
		 *  -1 = real beta^2
		 *   0 = complex beta^2, first of the complex-conjugate pair
		 *  +1 = complex beta^2, second of the complex-conjugate pair
		 *
		 * This info is needed in activate_mode().
		 */
		if(beta2.imag() == 0.0) {
		    modetypes.push_back(-1);
		} else {
		    modetypes.push_back(modetype++);
		    if(modetype == 2) {
			modetype = 0;
		    }
		}
	    }

	}

	size_t getNumEigenpairs () const final
        {
	    return num_eigenpairs;
	}

	std::complex<double> getEigenvalue (size_t i) const final
	{
	    return pEval[i];
	}

	std::vector<std::complex<double>> getEigenvector (size_t i) final
	{
	    if(is_matrix_real) {
		std::vector<std::complex<double>> result(n, 0.0);
		auto pBegin = pEvecReal.get() + i*n;
		get_complex_ht_sub(modetypes[i], &result[0], pBegin, n, n);
		return result;
	    } else {
		auto pBegin = pEvecComplex.get() + i*n;
		return std::vector<std::complex<double>>( pBegin, pBegin + n );
	    }
	}

    };

} // anonymous namespace


std::unique_ptr<ISolverResult> eigenSolver (
    std::unique_ptr< sparse_matrix<std::complex<double>> > matrix,
    bool is_matrix_real,
    double sigma,
    int num_eigenvalues
    )
{
    return make_unique<ArpackSolverResult>(std::move(matrix), is_matrix_real, sigma, num_eigenvalues);
}
