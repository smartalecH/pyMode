
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

#if defined(HAVE_MPI) && defined(HAVE_PETSC)
#include <mpi.h>
#include <petscsys.h>
#endif

#include <complex>
#include <list>
#include <stdexcept>
#include <cstring> // memset
#include <numeric> // accumulate

using std::complex;

#include <unistd.h>

#include <boost/throw_exception.hpp>

#include "wgms3d.h"
#include "mgp.h"
#include "stencil.h"
#include "sparse.h"
#include "pml.h"
#include "diffops.h"
#include "fortran_interface.h"
#include "solver.h"
#include "make_unique.h"
#include "simulation_parameters.h"

#define PML_POWER 2

namespace
{

    /// free-space impedance
    const double Z0 = 4e-7 * M_PI * 299792458;

    std::unique_ptr<sparse_matrix<complex<double>>>
    matrix_remove_zero_points (
	sparse_matrix<complex<double> > &in,
	const std::vector<int> &retain_list,
	unsigned number_of_unknowns
	)
    {
	auto out = std::unique_ptr<sparse_matrix<complex<double>>>( new sparse_matrix<complex<double>>(number_of_unknowns) );
	unsigned int inrow, j;

	in.order(1);
	for(inrow = 0; inrow < in.n; inrow++) {
	    int outrow = retain_list[inrow];
	    if(outrow >= 0) {
		for(j = in.indextable[inrow]; j < in.indextable[inrow+1]; j++) {
		    int incol = in.entries[j].j;
		    int outcol = retain_list[incol];
		    if(outcol >= 0) {
			out->add_entry(outrow, outcol, in.entries[j].v);
		    }
		}
	    }
	}

	return out;
    }

    /// Copy eigenvector, inserting known-zero points
    void unpackEigenvector (
	int num,
	const int *list,
	const std::complex<double> *src,
	std::complex<double> *dst
	)
    {
	int start, end = -1;
	while(end < num - 1) {
	    /* search start of next block to copy */
	    for(start = end + 1; start < num; start++) {
		if(list[start] >= 0) {
		    break;
		}
	    }
	    if(start == num) {
		/* there's no next block. */
		break;
	    }

	    /* search end of the block */
	    for(end = start + 1; end < num; end++) {
		if(list[end] < 0) {
		    break;
		}
	    }

	    /* copy data from start to end-1 */
	    int inc = 1;
	    // which would be faster? BLAS copy or std::copy?
	    COPY(end-start, (complex<double>*)src + list[start], inc,
		 dst + start, inc);
	}
    }
    
    /* The following function replaces references to points 'from' with
     * references to points 'to'. */
    template <class T>
    void handle_ghost_points_sub (T *m,
				  const int *from,
				  const int *to,
				  int Htan,
				  int Hnorm,
				  int sym,
				  int num)
    {
	double htan_symtype, hnorm_symtype;

	if(sym == 0) {
	    /* electric wall: Hnorm = 0 ("antisymmetric") */
	    htan_symtype = +1.0;
	    hnorm_symtype = -1.0;
	} else {
	    /* magnetic wall: Htan = 0 */
	    htan_symtype = -1.0;
	    hnorm_symtype = +1.0;
	}

	for(int k = 0; k < num; k++) {
	    int s_from = from[k];
	    int s_to   = to[k];
	    m[s_to + Htan]    += htan_symtype  * m[s_from + Htan];
	    m[s_to + Hnorm]   += hnorm_symtype * m[s_from + Hnorm];
	    m[s_from + Htan]   = 0.0;
	    m[s_from + Hnorm]  = 0.0;
	}
    }

    void
    grid_add_ghost_points (std::vector<double> &grid)
    {
	double delta;

	assert(grid.size() >= 2);

	delta = grid[1] - grid.front();
	grid.insert(grid.begin(), grid.front() - delta);

	delta = grid.back() - grid[grid.size()-2];
	grid.push_back(grid.back() + delta);
    }

    void init_pml_arrays_sub (const PML *pml,
			      const std::vector<double> &x,
			      int offset,
			      int nx,
			      /* store (1/s), (s'/s^3) and (\tilde\rho) here: */
			      std::vector<complex<double>> &s1,
			      std::vector<complex<double>> &s2,
			      std::vector<complex<double>> &sx)
    {
	int i;
	complex<double> s, sp;

	for(i = offset; i <= offset + nx; i++) {
	    s = pml->get_s(x[i]);
	    sp = pml->get_sprime(x[i]);
	    sx[i] = pml->get_stretched_x(x[i]);
	    s1[i] = 1.0 / s;
	    s2[i] = sp / (s * s * s);
	}
    }

    void AXPY2 (int n,
		complex<double> alpha,
		complex<double> *x1,
		complex<double> *x2,
		complex<double> *y)
    {
	while(n--) {
	    *y++ += alpha * (*x1++) * (*x2++);
	}
    }

} // anonymous namespace

namespace wgms3d {

    std::complex<double> Mode::getEffectiveIndex () const
    {
	return beta / sp->k0;
    }

    double Mode::getAlphaPer90Deg () const
    {
	if(sp->c == 0.0) {
	    return 0.0;
	} else {
	    return getAlphaPerUol() * M_PI / 2.0 * (1.0 / sp->c);
	}
    }

    ModeSolver::ModeSolver (void)
	: sp(std::make_shared<SimulationParameters>())
    {
    }

    ModeSolver::~ModeSolver (void)
    {
    }

    void ModeSolver::set_fd_mode (FDMode mode)
    {
	sp->fd_mode = mode;
    }

    FDMode ModeSolver::get_fd_mode (void) const
    {
	return sp->fd_mode;
    }

    bool ModeSolver::is_core_layer_defined (void) const
    {
	return sp->mgp->is_core_layer_defined();
    }

    void ModeSolver::set_wavelength (double lambda)
    {
	sp->k0 = 2.0 * M_PI / lambda;
    }

    void ModeSolver::set_curvature (double curvature)
    {
	sp->c = curvature;
    }

    void ModeSolver::set_walltype (int which_wall, int type)
    {
	assert(which_wall >= 0 && which_wall <= 3);
	assert(type == 0 || type == 1);

	sp->bconds[which_wall] = type;
    }

    void ModeSolver::set_five_point_standard (bool to)
    {
	sp->use_five_point_standard = to;
    }

    void
    ModeSolver::set_geometry (const char *mgp_filename)
    {
	sp->mgp = make_unique<MGP>(mgp_filename);
    }

/* -- GRID + GHOST-POINT HANDLING ---------------------------------------------------------- */

#if NUM_GHOST_POINTS == 1
    void
    ModeSolver::handle_ghost_points (complex<double> *m,
				     int ri,
				     int zi)
    {
	static const int Wpoints[3] = { 6, 7, 8 };
	static const int Vpoints[3] = { 5, 0, 1 };
	static const int Epoints[3] = { 4, 3, 2 };
	static const int Npoints[3] = { 8, 1, 2 };
	static const int Hpoints[3] = { 7, 0, 3 };
	static const int Spoints[3] = { 6, 5, 4 };

	int Htan, Hnorm;

	/* Vertical boundaries: tangential H is Hy, normal H is Hx */
	Htan = NSP; Hnorm = 0;
	if(ri == 1) {
	    handle_ghost_points_sub(m, Wpoints, sp->bcondsym[0] ? Vpoints : Epoints,
				    Htan, Hnorm, sp->bconds[0], 3);
	}
	if(ri == sp->nir) {
	    handle_ghost_points_sub(m, Epoints, sp->bcondsym[1] ? Vpoints : Wpoints,
				    Htan, Hnorm, sp->bconds[1], 3);
	}

	/* Horizontal boundaries: tangential H is Hx, normal H is Hy */
	Htan = 0; Hnorm = NSP;
	if(zi == 1) {
	    handle_ghost_points_sub(m, Spoints, sp->bcondsym[3] ? Hpoints : Npoints,
				    Htan, Hnorm, sp->bconds[3], 3);
	}
	if(zi == sp->niz) {
	    handle_ghost_points_sub(m, Npoints, sp->bcondsym[2] ? Hpoints : Spoints,
				    Htan, Hnorm, sp->bconds[2], 3);
	}
    }
#endif

    void
    ModeSolver::setGrid (const std::vector<double> & rho_grid,
			 const std::vector<double> & z_grid)
    {
	sp->_rhos = rho_grid;
	sp->_zs = z_grid;

	sp->nir = sp->_rhos.size();
	sp->niz = sp->_zs.size();
	sp->fcsize = sp->nir * sp->niz;

	grid_add_ghost_points(sp->_rhos);
	grid_add_ghost_points(sp->_zs);
    }

    const std::complex<double> *
    ModeSolver::get_stretched_rhos (void)
    {
	return &sp->stretched_rhos[1];
    }

    const std::complex<double> *
    ModeSolver::get_stretched_zs (void)
    {
	return &sp->stretched_zs[1];
    }

/* -- CODE FOR SETTING UP THE FINITE-DIFFERENCE SYSTEM MATRIX ---------------------------------- */

    void
    ModeSolver::get_epsis (const std::complex<double> **epsis,
			   int *ldepsis)
    {
	*epsis = this->epsis.get();
	*ldepsis = sp->_rhos.size();
    }

    void
    ModeSolver::add_matrix_entries (sparse_matrix<complex<double> > &A,
				    int to,      /* # of equation where we have to add the entries */
				    int Poffset, /* index of central stencil point */
				    const complex<double> *m)
    {
	auto nir = sp->nir;
	A.add_entry(to, Poffset,       *m++); /* P  */
	A.add_entry(to, Poffset+nir,   *m++); /* N  */
	A.add_entry(to, Poffset+nir+1, *m++); /* NE */
	A.add_entry(to, Poffset    +1, *m++); /* E  */
	A.add_entry(to, Poffset-nir+1, *m++); /* SE */
	A.add_entry(to, Poffset-nir,   *m++); /* S  */
	A.add_entry(to, Poffset-nir-1, *m++); /* SW */
	A.add_entry(to, Poffset    -1, *m++); /* W  */
	A.add_entry(to, Poffset+nir-1, *m); /* NW */
    }

    /* the following is just to store the command-line arguments. */
    struct pml_spec {
	char where;
	int numcells;
	double sigmamult;
    };

    void
    ModeSolver::add_pml (char where,
			 int numcells,
			 double sigmamult)
    {
	pml_spec newpml;
	newpml.where = where;
	newpml.numcells = numcells;
	newpml.sigmamult = sigmamult;
	pml_specs.push_back(newpml);
    }

    std::unique_ptr<sparse_matrix<complex<double>>>
    ModeSolver::initmatrix (Diffops & diffops,
			    double &nmax)
    {
	unsigned int i, j, k;

	/* ------------------- Setup PML configuration ----------------- */

	const auto & _rhos = sp->_rhos;
	const auto & _zs = sp->_zs;

	std::vector<complex<double>> s1(_rhos.size(), 1.0);
	std::vector<complex<double>> s2(_rhos.size(), 0.0);
	std::vector<complex<double>> t1(_zs.size(), 1.0);
	std::vector<complex<double>> t2(_zs.size(), 0.0);

	sp->stretched_rhos = std::vector<complex<double>>(_rhos.cbegin(), _rhos.cend());
	sp->stretched_zs = std::vector<complex<double>>(_zs.cbegin(), _zs.cend());

	PML *pml;
	for(pml_spec &p : pml_specs) {
	    switch(p.where) {
	    case 'n':
		pml = &sp->pml[0];
		pml->init(_zs[_zs.size()-1-p.numcells], 1, PML_POWER);
		pml->set_optimal_strength(_zs[_zs.size()-1-p.numcells+1] - _zs[_zs.size()-1-p.numcells],
					  _zs[_zs.size()-1] - _zs[_zs.size()-1-p.numcells],
					  sp->k0, p.sigmamult);
		init_pml_arrays_sub(pml, _zs, _zs.size()-1-p.numcells, p.numcells, t1, t2, sp->stretched_zs);
		break;
	    case 'e':
		pml = &sp->pml[1];
		pml->init(_rhos[_rhos.size()-1-p.numcells], 1, PML_POWER);
		pml->set_optimal_strength(_rhos[_rhos.size()-1-p.numcells+1] - _rhos[_rhos.size()-1-p.numcells],
					  _rhos[_rhos.size()-1] - _rhos[_rhos.size()-1-p.numcells],
					  sp->k0, p.sigmamult);
		init_pml_arrays_sub(pml, _rhos, _rhos.size()-1-p.numcells, p.numcells, s1, s2, sp->stretched_rhos);
		break;
	    case 's':
		pml = &sp->pml[2];
		pml->init(_zs[p.numcells], -1, PML_POWER);
		pml->set_optimal_strength(_zs[p.numcells] - _zs[p.numcells-1],
					  _zs[p.numcells] - _zs[0],
					  sp->k0, p.sigmamult);
		init_pml_arrays_sub(pml, _zs, 0, p.numcells, t1, t2, sp->stretched_zs);
		break;
	    case 'w':
		pml = &sp->pml[3];
		pml->init(_rhos[p.numcells], -1, PML_POWER);
		pml->set_optimal_strength(_rhos[p.numcells] - _rhos[p.numcells-1],
					  _rhos[p.numcells] - _rhos[0],
					  sp->k0, p.sigmamult);
		init_pml_arrays_sub(pml, _rhos, 0, p.numcells, s1, s2, sp->stretched_rhos);
		break;
	    }
	}

	/* -------------------------------------------------------------------------- */

	epsis = sp->mgp->get_epsis_on_grid(_rhos, _zs);
	const unsigned ldepsis = _rhos.size();

	nmax = 1.0;

	std::cout << "Setting up FD system matrix (initial dimension = "
		  << 2*sp->fcsize << ")... " << std::endl;

	auto A0 = make_unique<sparse_matrix<complex<double>>>(2*sp->fcsize);

	const auto nir = sp->nir;
	for(j = NUM_GHOST_POINTS; j < _zs.size() - NUM_GHOST_POINTS; j++) {
	    const double zp = _zs[j];
	    const double n = _zs[j+1]-zp;
	    const double s = zp-_zs[j-1];

	    for(i = NUM_GHOST_POINTS; i < _rhos.size() - NUM_GHOST_POINTS; i++) {
		const double rp = _rhos[i];
		const double e = _rhos[i+1]-rp;
		const double w = rp-_rhos[i-1];

		/* number of current grid point */
		/* The current grid point has the number gpn = '(j-1)*nir
		 * + (i-1)'. The FD equation for the Hrho field component is
		 * stored in matrix row gpn, while the equation for the Hz
		 * field is stored in matrix row gpn+fcsize. */
		const int gpn = (j-NUM_GHOST_POINTS)*nir + (i-NUM_GHOST_POINTS);

		complex<double> *M0;

		complex<double> epsp = epsis[j*ldepsis + i];
		complex<double> complex_n = sqrt(epsp);
		if(complex_n.real() > nmax) {
		    nmax = complex_n.real();
		}

		if(debugmgp) {
		    std::cout << "* (" << i << "," << j << "): n = "
			      << complex_n << std::endl;
		}

		bool standard;

		const direction dirs[] = DIRS;
		M0 = diffops.calculate_diffop(rp, zp, epsp, dirs);
		if(!M0) {
		    /* No interfaces found in this stencil. Use explicit
		     * expressions for FD weights => less "numerical
		     * noise", less non-zero matrix entries. */
		    M0 = diffops.get_standard_diffop(n, e, s, w);
		    standard = true;
		} else {
		    standard = false;
		}

		if(s1[i] != 1.0 || t1[j] != 1.0) {
		    /* We're inside a PML. Replace the derivatives by the
		     * "PML-stretched" versions, so that in the remainder
		     * of the program we can simply pretend there is no
		     * PML (using the complex stretched \rho in the
		     * differential equations, though) -- everything else
		     * is absorbed in these matrices. */
		    if(standard) {
			/* the M0 returned by get_standard_diffop() must
			 * not be modified, so we make a backup here. */
			complex<double> *M0orig = M0;
			M0 = new complex<double>[(2*NDO) * (2*NSP)];
			memcpy(M0, M0orig, (2*NDO) * (2*NSP) * sizeof(*M0));
			standard = false;
			/* (The memory just allocated for M0 will be freed
			 * in the Diffops destructor.) */
		    }
		    for(k = 0; k < 2; k++) {
			/* h^X_\rho\rho  ==>  (1/s^2) h^X_\rho\rho - (s'/s^3) h^X_\rho */
			SCAL(2*NSP, s1[i]*s1[i], M0 + 2 + k*NDO, 2*NDO);
			AXPY(2*NSP, -s2[i],      M0 + 0 + k*NDO, 2*NDO,
			     M0 + 2 + k*NDO, 2*NDO);
			/* h^X_zz        ==>  (1/t^2) h^X_zz - (t'/t^3) h^X_z */
			SCAL(2*NSP, t1[j]*t1[j], M0 + 4 + k*NDO, 2*NDO);
			AXPY(2*NSP, -t2[j],      M0 + 1 + k*NDO, 2*NDO,
			     M0 + 4 + k*NDO, 2*NDO);
			/* h^X_\rho      ==>  (1/s) h^X_\rho */
			SCAL(2*NSP, s1[i],       M0 + 0 + k*NDO, 2*NDO);
			/* h^X_z         ==>  (1/t) h^X_z */
			SCAL(2*NSP, t1[j],       M0 + 1 + k*NDO, 2*NDO);
			/* h^X_\rho z    ==>  (1/(s*t)) h^X_\rho z */
			SCAL(2*NSP, s1[i]*t1[j], M0 + 3 + k*NDO, 2*NDO);
		    }
		}

		/* If not standard, save the derivatives for later use in
		 * export_mode(). */
		if(!standard) {
		    diffops.store_diffops(M0, i, j);
		}

		/* Now continue to setup system matrix. */
		const complex<double> onepxc = 1.0 + sp->stretched_rhos[i] * sp->c;
		const complex<double> onepxc2 = onepxc * onepxc;
		complex<double> coeffs1[2*NSP];
		complex<double> coeffs2[2*NSP];
		std::memset(coeffs1, 0, sizeof(coeffs1));
		std::memset(coeffs2, 0, sizeof(coeffs2));

		/* Scalar Helmholtz operator: */
		AXPY(2*NSP, onepxc2, M0 + 2, 2*NDO, coeffs1, 1);
		AXPY(2*NSP, onepxc2, M0 + 4, 2*NDO, coeffs1, 1);
		coeffs1[0]   += onepxc2*sp->k0*sp->k0*epsp;

		AXPY(2*NSP, onepxc2, M0 + NDO+2, 2*NDO, coeffs2, 1);
		AXPY(2*NSP, onepxc2, M0 + NDO+4, 2*NDO, coeffs2, 1);
		coeffs2[NSP] += onepxc2*sp->k0*sp->k0*epsp;

		if(sp->c != 0.0) {
		    /* --- Additional curvature terms --- */
		    AXPY(2*NSP, 3*sp->c*onepxc, M0 + 0,     2*NDO, coeffs1, 1);
		    AXPY(2*NSP, 2*sp->c*onepxc, M0 + NDO+1, 2*NDO, coeffs1, 1);
		    coeffs1[0] += sp->c*sp->c;
		    AXPY(2*NSP, sp->c*onepxc,   M0 + NDO+0, 2*NDO, coeffs2, 1);
		}

		/* Make sure that references to ghost points are replaced
		 * by references to real grid points (either inner grid
		 * points or boundary points). */
		handle_ghost_points(coeffs1, i, j);
		handle_ghost_points(coeffs2, i, j);

		add_matrix_entries(*A0, gpn, gpn, coeffs1);
		add_matrix_entries(*A0, gpn, gpn+sp->fcsize, coeffs1 + NSP);
		add_matrix_entries(*A0, gpn+sp->fcsize, gpn, coeffs2);
		add_matrix_entries(*A0, gpn+sp->fcsize, gpn+sp->fcsize, coeffs2 + NSP);
	    }
	}

	std::cout << "Stored "
		  << diffops.get_num_stored_diffops() << "/" << sp->fcsize
		  << " non-standard diffops (~"
		  << (diffops.get_num_stored_diffops()*(2*NSP)*(2*NDO)*sizeof(complex<double>))/(1<<20)
		  << "MB)." << std::endl;

	return A0;
    }

#if 0
    void
    wgms3d_mode::adjust_phase (void)
    {
	int i, j;
	double avgreal = 0.0;
	double avgimag = 0.0;
	int n = 0;

	wg->activate_mode(this);
    }
#endif

    /// One instance of this class will be produced by wgms3d and then
    /// passed in a shared_ptr to the the various
    /// VectorModes. (current design: we assume that it is known
    /// before the solving process which of the derived fields are to
    /// be calculated [= third argument to calculateModes()]; only the
    /// corresponding matrices and vectors will then be set up by
    /// ModeSolver.).
    class DerivedFieldMatrices
    {
    public:
	unsigned which;
	sparse_matrix<std::complex<double>> matrix_Hr_to_Er;
	sparse_matrix<std::complex<double>> matrix_Hz_to_Er;
	sparse_matrix<std::complex<double>> matrix_Hr_to_Ez;
	sparse_matrix<std::complex<double>> matrix_Hz_to_Ez;
	sparse_matrix<std::complex<double>> matrix_Hr_to_Ep;
	sparse_matrix<std::complex<double>> matrix_Hz_to_Ep;
	sparse_matrix<std::complex<double>> matrix_Hr_to_Hp;
	sparse_matrix<std::complex<double>> matrix_Hz_to_Hp;
	std::unique_ptr< std::complex<double>[] > vector_Hz_to_Er;
	std::unique_ptr< std::complex<double>[] > vector_Hr_to_Ez;
    };

    // TODO maybe move somewhere else?
    void ModeSolver::initDerivedFieldMatrices (Diffops & diffops,
					       unsigned which)
    {
	auto newdfm = std::make_shared<DerivedFieldMatrices>();
	deriv_field_mats = newdfm;

	auto &dfm = *newdfm;
	dfm.which = which;

	if(!which) {
	    return;
	}

	static const complex<double> jay(0,1);


	if(which & 1) {
	    dfm.matrix_Hr_to_Er.init(sp->fcsize, sp->fcsize);
	    dfm.matrix_Hz_to_Er.init(sp->fcsize, sp->fcsize);
	    dfm.vector_Hz_to_Er.reset(new complex<double>[sp->fcsize]);
	}

	if(which & 2) {
	    dfm.matrix_Hr_to_Ez.init(sp->fcsize, sp->fcsize);
	    dfm.matrix_Hz_to_Ez.init(sp->fcsize, sp->fcsize);
	    dfm.vector_Hr_to_Ez.reset(new complex<double>[sp->fcsize]);
	}

	if(which & 4) {
	    dfm.matrix_Hr_to_Ep.init(sp->fcsize, sp->fcsize);
	    dfm.matrix_Hz_to_Ep.init(sp->fcsize, sp->fcsize);
	}

	if(which & 8) {
	    dfm.matrix_Hr_to_Hp.init(sp->fcsize, sp->fcsize);
	    dfm.matrix_Hz_to_Hp.init(sp->fcsize, sp->fcsize);
	}

	const auto & _rhos = sp->_rhos;
	const auto & _zs = sp->_zs;
	const unsigned ldepsis = _rhos.size();

	for(unsigned j = NUM_GHOST_POINTS; j < _zs.size() - NUM_GHOST_POINTS; j++) {
	    for(unsigned i = NUM_GHOST_POINTS; i < _rhos.size() - NUM_GHOST_POINTS; i++) {
		const int gpn = (j-NUM_GHOST_POINTS)*sp->nir + (i-NUM_GHOST_POINTS);
		complex<double> eps = epsis[j*ldepsis + i];
		complex<double> *M0 = diffops.get_diffops(_rhos, _zs, i, j);
		const complex<double> onepxc = 1.0 + sp->stretched_rhos[i]*sp->c;
		const complex<double> Zkn2 = Z0 / (sp->k0 * eps);
		complex<double> scale;
		complex<double> coeffs[2*NSP];

		if(which & 1) {
		    dfm.vector_Hz_to_Er[gpn] = -Zkn2 / onepxc; /* h^z */
		    std::memset(coeffs, 0, sizeof(coeffs));
		    scale = +Zkn2 * onepxc;
		    AXPY(2*NSP, scale, M0 + 3,       2*NDO, coeffs, 1); /* h^\rho_{\rho z} */
		    AXPY(2*NSP, scale, M0 + NDO + 4, 2*NDO, coeffs, 1); /* h^z_{z z} */
		    scale = +Zkn2 * sp->c;
		    AXPY(2*NSP, scale, M0 + 1,       2*NDO, coeffs, 1); /* h^\rho_z */
		    handle_ghost_points(coeffs, i, j);
		    add_matrix_entries(dfm.matrix_Hr_to_Er, gpn, gpn, coeffs + 0);
		    add_matrix_entries(dfm.matrix_Hz_to_Er, gpn, gpn, coeffs + NSP);
		}

		if(which & 2) {
		    dfm.vector_Hr_to_Ez[gpn] = Zkn2 / onepxc; /* h^\rho */
		    std::memset(coeffs, 0, sizeof(coeffs));
		    coeffs[0] = -Zkn2*sp->c*sp->c / onepxc; /* h^\rho */
		    scale = -Zkn2 * onepxc;
		    AXPY(2*NSP, scale, M0 + 2,       2*NDO, coeffs, 1); /* h^\rho_{\rho\rho} */
		    AXPY(2*NSP, scale, M0 + NDO + 3, 2*NDO, coeffs, 1); /* h^z_{\rho z} */
		    scale = -3.0 * sp->c * Zkn2;
		    AXPY(2*NSP, scale, M0 + 0,       2*NDO, coeffs, 1); /* h^\rho_\rho */
		    scale = -2.0 * sp->c * Zkn2;
		    AXPY(2*NSP, scale, M0 + NDO + 1, 2*NDO, coeffs, 1); /* h^z_z */
		    handle_ghost_points(coeffs, i, j);
		    add_matrix_entries(dfm.matrix_Hr_to_Ez, gpn, gpn, coeffs + 0);
		    add_matrix_entries(dfm.matrix_Hz_to_Ez, gpn, gpn, coeffs + NSP);
		}

		if(which & 4) {
		    std::memset(coeffs, 0, sizeof(coeffs));
		    scale = jay * Zkn2;
		    AXPY(2*NSP, scale, M0 + 1,       2*NDO, coeffs, 1); /* h^\rho_z */
		    scale *= -1.0;
		    AXPY(2*NSP, scale, M0 + NDO + 0, 2*NDO, coeffs, 1); /* h^z_\rho */
		    handle_ghost_points(coeffs, i, j);
		    add_matrix_entries(dfm.matrix_Hr_to_Ep, gpn, gpn, coeffs + 0);
		    add_matrix_entries(dfm.matrix_Hz_to_Ep, gpn, gpn, coeffs + NSP);
		}

		if(which & 8) {
		    std::memset(coeffs, 0, sizeof(coeffs));
		    coeffs[0] = jay * sp->c; /* h^\rho */
		    scale = jay * onepxc;
		    AXPY(2*NSP, scale, M0 + 0,       2*NDO, coeffs, 1); /* h^\rho_\rho */
		    AXPY(2*NSP, scale, M0 + NDO + 1, 2*NDO, coeffs, 1); /* h^z_z */
		    handle_ghost_points(coeffs, i, j);
		    add_matrix_entries(dfm.matrix_Hr_to_Hp, gpn, gpn, coeffs + 0);
		    add_matrix_entries(dfm.matrix_Hz_to_Hp, gpn, gpn, coeffs + NSP);
		}
	    }
	}

    }

    void ModeSolver::calculateModes (
	unsigned number_of_modes,
	double near_this_effective_index,
	unsigned derived_field_mask
	)
    {
	int rank = 0;
#if defined(HAVE_MPI)
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif

	std::unique_ptr<sparse_matrix<complex<double>>> A;
	double shift_invert_sigma;

	Diffops diffops(sp->mgp.get(), sp);

	if(rank == 0)
	{
	    // The root process prepares the FD matrix.
	    double nmax;

	    /* Prepare the finite-difference system matrix. */
	    A = initmatrix(diffops, nmax);

	    /* Remove those unknowns from the system. (We know the
	     * field must be zero there, so this adds unnecessary
	     * information to the system, which can result in spurious
	     * eigenvalues.) - TODO: do this in one step? */
	    A = matrix_remove_zero_points(*A, sp->getRetainList(), sp->getNumUnknowns());

	    std::cout << "Final matrix dimension is "
		      << A->n << "; "
		      << A->length << " non-zero entries." << std::endl;

#if defined(HAVE_MPI)
	    {
		int n = A->n;
		MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	    }
#endif

	    /* Prepare for ARPACK's shift-and-invert mode. */
	    if(near_this_effective_index < 0.0) {
		near_this_effective_index = nmax;
	    }
	    std::cout << "Searching for modes near n_eff = "
		      << near_this_effective_index << "." << std::endl;
	    shift_invert_sigma = std::pow(sp->k0*near_this_effective_index, 2);
#if defined(HAVE_MPI)
	    MPI_Bcast(&shift_invert_sigma, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
#endif
	}
	else
	{
#if defined(HAVE_MPI)
	    // The non-root processes only need to know the matrix
	    // dimension and the sigma.
	    int n;
	    MPI_Bcast(&n, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	    MPI_Bcast(&shift_invert_sigma, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	    A.reset( new sparse_matrix<complex<double>>(n) );
#else
	    // shouldn't be reached. suppress clang warning:
	    shift_invert_sigma = -1.0;	    
#endif
	}

	//std::cout << "A->n = " << A->n << ", A->length = " << A->length << ", sigma = " << shift_invert_sigma << std::endl;
	const bool real_matrix = !sp->mgp->hasComplexIndices() && pml_specs.size() == 0;
	solver_result = eigenSolver(std::move(A), real_matrix, shift_invert_sigma, number_of_modes);

	if(rank == 0) {
	    // We initialize these matrices at this point (after
	    // eigensolution) to minimize the maximum memory usage.
	    initDerivedFieldMatrices(diffops, derived_field_mask);

	}
    }

    std::unique_ptr<Mode> ModeSolver::getMode (unsigned i)
    {
	//std::cout << "getMode(" << i << ")" << std::endl;
	auto beta2 = solver_result->getEigenvalue(i); /* This is beta^2 */
	auto migamma = sqrt(beta2); /* sqrt(beta^2) = (-i * gamma) */

	/* Make sure we have modes with non-negative real part of the
	 * effective index = "forward-travelling waves". */
	if(migamma.real() < 0) {
	    migamma = -migamma;
	}
    
	if(sp->fd_mode == FDMode::Scalar)
	{
	    return make_unique<ScalarMode>(i, migamma, sp, solver_result);
	}
	else
	{
	    return make_unique<VectorMode>(i, migamma, sp, solver_result, deriv_field_mats);
	}
    }

    Mode::Mode (unsigned number,
		std::complex<double> beta,
		std::shared_ptr<SimulationParameters> simulation_parameters)
	: number(number), beta(beta), sp(simulation_parameters)
    {
	rank = 0;
#if defined(HAVE_MPI)
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif
    }

    char VectorMode::estimatePolarization () const
    {
	ensurePrimaryFields();

	if(rank != 0) {
	    return '-';
	}

	auto sumAbs = [](double partial_sum, const std::complex<double> & val) {
	    return partial_sum + std::fabs(val);
	};
	auto sum_Hr = std::accumulate(Hr.begin(), Hr.end(), double(0.0), sumAbs);
	auto sum_Hz = std::accumulate(Hz.begin(), Hz.end(), double(0.0), sumAbs);

	if(sum_Hr > 2*sum_Hz) {
	    return 'V';
	} else if(sum_Hz > 2*sum_Hr) {
	    return 'H';
	} else {
	    return '?';
	}
    }

    const std::vector<std::complex<double>> & VectorMode::getHr () const
    {
	ensurePrimaryFields();
	return Hr;
    }

    const std::vector<std::complex<double>> & VectorMode::getHz () const
    {
	ensurePrimaryFields();
	return Hz;
    }

    void VectorMode::ensurePrimaryFields () const
    {
	std::lock_guard<std::mutex> lock(mutables);

	if(Ht_valid) {
	    return;
	}

	//std::cout << rank << " Getting EV #" << number << std::endl;
	const auto eigenvector = solver_result->getEigenvector(number);

	if(rank == 0) {
	    const std::vector<int> & retain_list = sp->getRetainList();

	    Hr.resize(sp->fcsize, 0.0);
	    Hz.resize(sp->fcsize, 0.0);
	    unpackEigenvector(sp->fcsize, retain_list.data(), eigenvector.data(), Hr.data());
	    unpackEigenvector(sp->fcsize, retain_list.data()+sp->fcsize, eigenvector.data(), Hz.data());

	    /* Add phase factor to the entire mode field such that its
	     * phase averaged (in some sense) over the core region
	     * (defined by the 'C' directive in the geometry file) is
	     * zero. This does not have any physical meaning and is
	     * only for convenience when visualizing the mode
	     * fields. */
	    if(sp->mgp->is_core_layer_defined())
	    {
		// TODO optimize & make optional
		int ldHt = sp->nir;
		double avgreal = 0.0, avgimag = 0.0;
		unsigned n = 0;
		for(int j = 0; j < sp->niz; j++) {
		    double zz = sp->_zs[j+NUM_GHOST_POINTS];
		    for(int i = 0; i < sp->nir; i++) {
			double rr = sp->_rhos[i+NUM_GHOST_POINTS];
			if(sp->mgp->is_in_core(rr,zz)) {
			    avgreal += Hr[j*ldHt + i].real() + Hz[j*ldHt + i].real();
			    avgimag += Hr[j*ldHt + i].imag() + Hz[j*ldHt + i].imag();
			    n++;
			}
		    }
		}

		if(n != 0) {
		    double phase = atan2(avgimag / n, avgreal / n);
		    const std::complex<double> factor = exp(std::complex<double>(0.0, -phase));
		    SCAL(sp->fcsize, factor, Hr.data(), 1);
		    SCAL(sp->fcsize, factor, Hz.data(), 1);
		}
	    }

	}

	Ht_valid = true;
    }

    std::vector<std::complex<double>> ScalarMode::getField () const
    {
	const auto eigenvector = solver_result->getEigenvector(number);
	std::vector<std::complex<double>> field;

	if(rank == 0) {
	    const std::vector<int> & retain_list = sp->getRetainList();

	    field.resize(sp->fcsize, 0.0);
	    unpackEigenvector(sp->fcsize, retain_list.data()+sp->fcsize, eigenvector.data(), field.data());

#if 0
	  TODO:
	    /* Add phase factor to the entire mode field such that its
	     * phase averaged (in some sense) over the core region
	     * (defined by the 'C' directive in the geometry file) is
	     * zero. This does not have any physical meaning and is
	     * only for convenience when visualizing the mode
	     * fields. */
#endif

	}

	return field;
    }

    std::vector<std::complex<double>> VectorMode::computeDerivedField (
	unsigned which,
	std::function<void(std::complex<double>*)> f
	) const
    {
	// perform collective call:
	ensurePrimaryFields();

	std::vector<std::complex<double>> result;

	if(rank == 0) {
	    // only really do and return something if we're the MPI
	    // root process...
	    if(!(dfm->which & which)) {
		BOOST_THROW_EXCEPTION(std::runtime_error("Requested derived field that was not announced during calculateModes()."));
	    }
	    result.resize(sp->fcsize);
	    f(result.data());
	}

	return result;
    }

    std::vector<std::complex<double>> VectorMode::getEr () const
    {
	return computeDerivedField(1, [this](std::complex<double> *result) {
		dfm->matrix_Hr_to_Er.vecmult(result, Hr.data());
		dfm->matrix_Hz_to_Er.vecmultadd(result, Hz.data());
		SCAL(sp->fcsize, 1.0/beta, result, 1);
		AXPY2(sp->fcsize, beta, dfm->vector_Hz_to_Er.get(), Hz.data(), result);
	    });
    }

    std::vector<std::complex<double>> VectorMode::getEz () const
    {
	return computeDerivedField(2, [this](std::complex<double> *result) {
		dfm->matrix_Hr_to_Ez.vecmult(result, Hr.data());
		dfm->matrix_Hz_to_Ez.vecmultadd(result, Hz.data());
		SCAL(sp->fcsize, 1.0/beta, result, 1);
		AXPY2(sp->fcsize, beta, dfm->vector_Hr_to_Ez.get(), Hr.data(), result);
	    });
    }

    std::vector<std::complex<double>> VectorMode::getEp () const
    {
	return computeDerivedField(4, [this](std::complex<double> *result) {
		dfm->matrix_Hr_to_Ep.vecmult(result, Hr.data());
		dfm->matrix_Hz_to_Ep.vecmultadd(result, Hz.data());
	    });
    }

    std::vector<std::complex<double>> VectorMode::getHp () const
    {
	return computeDerivedField(8, [this](std::complex<double> *result) {
		dfm->matrix_Hr_to_Hp.vecmult(result, Hr.data());
		dfm->matrix_Hz_to_Hp.vecmultadd(result, Hz.data());
		SCAL(sp->fcsize, 1.0/beta, result, 1);
	    });
    }

    const std::vector<int> & SimulationParameters::getRetainList () const
    {
	ensureRetainList();
	return retain_list;
    }

    unsigned SimulationParameters::getNumUnknowns () const
    {
	ensureRetainList();
	return number_of_unknowns;
    }

    void SimulationParameters::ensureRetainList () const
    {
	if(retain_list.size() != 0) {
	    return;
	}

	std::lock_guard<std::mutex> lock(mutables);

	assert(fcsize > 0);
	retain_list.resize(2*fcsize, 0);

	auto prepare_retain_list_sub = [this](int ri, int rinc, int zi, int zinc, int n, int x_or_y) {
	    int offset = ri + zi*nir;
	    if(x_or_y == 1) {
		offset += fcsize;
	    }
	    while(n--) {
		retain_list[offset] = -1;
		offset += rinc + zinc*nir;
	    }
	};

	/* Eliminate boundary points with Dirichlet BCs. */
	if(bcondsym[0] == 0) {
	    /* W wall: Hnorm = Hx, Htan = Hy. */
	    prepare_retain_list_sub(0, 0, 0, 1, niz, bconds[0] == 1);
	}
	if(bcondsym[1] == 0) {
	    /* E wall: Hnorm = Hx, Htan = Hy. */
	    prepare_retain_list_sub(nir-1, 0, 0, 1, niz, bconds[1] == 1);
	}

	if(bcondsym[2] == 0) {
	    /* N wall: Hnorm = Hy, Htan = Hx. */
	    prepare_retain_list_sub(0, 1, niz-1, 0, nir, bconds[2] == 0);
	}
	if(bcondsym[3] == 0) {
	    /* S wall: Hnorm = Hy, Htan = Hx. */
	    prepare_retain_list_sub(0, 1, 0, 0, nir, bconds[3] == 0);
	}

	/* Eliminate points for semi-vectorial calculation. */
	int i;
	if(fd_mode == FDMode::SemiVectorialHz || fd_mode == FDMode::Scalar) {
	    /* only Hz field retained. */
	    for(i = 0; i < fcsize; i++) {
		retain_list[i] = -1;
	    }
	} else if(fd_mode == FDMode::SemiVectorialHr) {
	    /* only Hrho field retained. */
	    for(i = 0; i < fcsize; i++) {
		retain_list[fcsize + i] = -1;
	    }
	}

	/* Now give the retained points new numbers. */
	for(i = 0, number_of_unknowns = 0; i < 2*fcsize; i++) {
	    if(retain_list[i] == 0) {
		retain_list[i] = number_of_unknowns++;
	    }
	}

	std::cout << "Eliminated "
		  << 2*fcsize - number_of_unknowns
		  << " unknowns with Dirichlet BCs." << std::endl;

    }

} // namespace wgms3d
