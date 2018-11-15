
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2013  Michael Krause <m.krause@tu-harburg.de>

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

#include <iostream>
#include <list>
#include <cassert>
#include <cstring>
#include <cmath>
#include <complex>
#include <algorithm> // for_each

using std::complex;

#include "diffops.h"
#include "wgms3d.h"
#include "mgp.h"
#include "stencil.h"
#include "fortran_interface.h"
#include "simulation_parameters.h"

namespace
{

    /* Add entries to interface-equation matrix that are proportional
     * to h_1^- */
    void diffops_add_h1_entries (complex<double> *row,
				 const complex<double> h1coefficient)
    {
	assert(NDO == 5);

	/*  Coefficient for RmZ */
	row[2*12] +=  h1coefficient;
	/*  Coefficient for ZmR */
	row[7*12] += -h1coefficient;
    }

    /* Add entries to interface-equation matrix that are proportional
     * to h_2^- */
    void diffops_add_h2_entries (complex<double> *row,
				 const double C,
				 const double S,
				 const complex<double> h2coefficient)
    {
	assert(NDO == 5);

	/*  Coefficient for RmZZ */
	row[ 5*12] +=  C*h2coefficient;
	/*  Coefficient for ZmRZ */
	row[10*12] += -C*h2coefficient;
	/*  Coefficient for RmRZ */
	row[ 4*12] += -S*h2coefficient;
	/*  Coefficient for ZmRR */
	row[ 9*12] +=  S*h2coefficient;
    }

    void make_taylor_matrix (complex<double> *M,
			     double dr,
			     double dz)
    {
	int i;

	std::memset(M, 0, NV*NV*sizeof(M[0]));

	for(i = 0; i < NV; i += NF) { /* This is the same for the R and Z field components */
	    M[0+i + (0+i)*NV] = 1.0;
	    M[0+i + (1+i)*NV] = dr;
	    M[0+i + (2+i)*NV] = dz;
	    M[0+i + (3+i)*NV] = (dr*dr)/2.0;
	    M[0+i + (4+i)*NV] = dr*dz;
	    M[0+i + (5+i)*NV] = (dz*dz)/2.0;

	    M[1+i + (1+i)*NV] = 1.0;
	    M[1+i + (3+i)*NV] = dr;
	    M[1+i + (4+i)*NV] = dz;

	    M[2+i + (2+i)*NV] = 1.0;
	    M[2+i + (4+i)*NV] = dr;
	    M[2+i + (5+i)*NV] = dz;

	    M[3+i + (3+i)*NV] = 1.0;

	    M[4+i + (4+i)*NV] = 1.0;

	    M[5+i + (5+i)*NV] = 1.0;

#if NDO >= 9
	    M[0+i + (6+i)*NV] = (dr*dr*dr)/6.0;
	    M[0+i + (7+i)*NV] = (dr*dr*dz)/2.0;
	    M[0+i + (8+i)*NV] = (dr*dz*dz)/2.0;
	    M[0+i + (9+i)*NV] = (dz*dz*dz)/6.0;

	    M[1+i + (6+i)*NV] = (dr*dr)/2.0;
	    M[1+i + (7+i)*NV] = (dr*dz);
	    M[1+i + (8+i)*NV] = (dz*dz)/2.0;

	    M[2+i + (7+i)*NV] = (dr*dr)/2.0;
	    M[2+i + (8+i)*NV] = (dr*dz);
	    M[2+i + (9+i)*NV] = (dz*dz)/2.0;

	    M[3+i + (6+i)*NV] = dr;
	    M[3+i + (7+i)*NV] = dz;

	    M[4+i + (7+i)*NV] = dr;
	    M[4+i + (8+i)*NV] = dz;

	    M[5+i + (8+i)*NV] = dr;
	    M[5+i + (9+i)*NV] = dz;

	    M[6+i + (6+i)*NV] = 1.0;

	    M[7+i + (7+i)*NV] = 1.0;

	    M[8+i + (8+i)*NV] = 1.0;

	    M[9+i + (9+i)*NV] = 1.0;
#endif

#if NDO >= 14
	    M[0+i + (10+i)*NV] = (dr*dr*dr*dr)/24.0;
	    M[0+i + (11+i)*NV] = (dr*dr*dr*dz)/6.0;
	    M[0+i + (12+i)*NV] = (dr*dr*dz*dz)/4.0;
	    M[0+i + (13+i)*NV] = (dr*dz*dz*dz)/6.0;
	    M[0+i + (14+i)*NV] = (dz*dz*dz*dz)/24.0;

	    M[1+i + (10+i)*NV] = (dr*dr*dr)/6.0;
	    M[1+i + (11+i)*NV] = (dr*dr*dz)/2.0;
	    M[1+i + (12+i)*NV] = (dr*dz*dz)/2.0;
	    M[1+i + (13+i)*NV] = (dz*dz*dz)/6.0;

	    M[2+i + (11+i)*NV] = (dr*dr*dr)/6.0;
	    M[2+i + (12+i)*NV] = (dr*dr*dz)/2.0;
	    M[2+i + (13+i)*NV] = (dr*dz*dz)/2.0;
	    M[2+i + (14+i)*NV] = (dz*dz*dz)/6.0;

	    M[3+i + (10+i)*NV] = (dr*dr)/2.0;
	    M[3+i + (11+i)*NV] = (dr*dz);
	    M[3+i + (12+i)*NV] = (dz*dz)/2.0;

	    M[4+i + (11+i)*NV] = (dr*dr)/2.0;
	    M[4+i + (12+i)*NV] = (dr*dz);
	    M[4+i + (13+i)*NV] = (dz*dz)/2.0;

	    M[5+i + (12+i)*NV] = (dr*dr)/2.0;
	    M[5+i + (13+i)*NV] = (dr*dz);
	    M[5+i + (14+i)*NV] = (dz*dz)/2.0;

	    M[6+i + (10+i)*NV] = dr;
	    M[6+i + (11+i)*NV] = dz;

	    M[7+i + (11+i)*NV] = dr;
	    M[7+i + (12+i)*NV] = dz;

	    M[8+i + (12+i)*NV] = dr;
	    M[8+i + (13+i)*NV] = dz;

	    M[9+i + (13+i)*NV] = dr;
	    M[9+i + (14+i)*NV] = dz;

	    M[10+i + (10+i)*NV] = 1.0;

	    M[11+i + (11+i)*NV] = 1.0;

	    M[12+i + (12+i)*NV] = 1.0;

	    M[13+i + (13+i)*NV] = 1.0;

	    M[14+i + (14+i)*NV] = 1.0;
#endif
	}
    }

    void delete_diffop (std::pair<const int, std::complex<double>*> &pair)
    {
	delete[] pair.second;
    }

} // anonymous namespace

namespace wgms3d {

    void
    Diffops::make_curv_interface_matrix (
	complex<double> *MLR,
	double theta,
	double d,
	complex<double> m,
	complex<double> p,
	double rho, double z)
    {
	const double C = cos(theta);
	const double S = sin(theta);
	const complex<double> Deps = p - m;
	/* The next line is a very simple way to implement a scalar
	 * computation. This is not in the paper. Might not be
	 * obvious. Derivations are in my notes from 20.5.2011. */
	const complex<double> D = (sp->fd_mode == FDMode::Scalar) ? 0 : Deps / m;
	int i;

#if 0
#warning "d set to zero"
	d = 0;
#endif

#if NDO != 5
#error "make_curv_interface_matrix doesn't support NDO>5"
#endif

	std::memset(MLR, 0, NV*NV*sizeof(MLR[0]));

	/* Initialize diagonal entries to unity (= trivial interface
	 * equations = no discontinuity) */
	for(i = 0; i < 12; i++) {
	    MLR[i+i*12] = 1.0;
	}

	/* Let's see if we're inside a PML region. If yes, use special
	 * interface code. */

	/* If we're in a z-directed PML: interfaces in such PMLs are
	 * not supported right now! */
	if(sp->pml[0].is_inside(z) || sp->pml[2].is_inside(z)) {
	    /* Equations (39)-(42) would have to be implemented
	     * here. I am too lazy right now. */
	    std::cerr << "No interfaces allowed in North or South PML regions." << std::endl;
	    exit(1);
	}

	/* Check if we are in a horizontal PML */
	PML *pml = NULL;
	if(sp->pml[1].is_inside(rho)) {
	    pml = sp->pml + 1;
	} else if(sp->pml[3].is_inside(rho)) {
	    pml = sp->pml + 3;
	}
	if(pml != NULL) {
	    /* We are inside a rho-directed PML. Allow only rho-oriented
	     * (horizontal) interfaces, and use appropriate boundary
	     * conditions. */
	    if(std::fabs(C) > 1e-10) {
		std::cerr << "Only horizontal interfaces allowed in East or West PML regions." << std::endl;
		exit(1);
	    }

	    /* Everything's fine. */
	    complex<double> s = pml->get_s(rho);
	    complex<double> ss = pml->get_sprime(rho) / (s*s);

	    /* RpZ */
	    MLR[2+2*12] += D;
	    MLR[2+7*12]  = -D/s;

	    /* RpRZ */
	    MLR[4+4*12] += D;
	    MLR[4+7*12]  = ss*D;
	    MLR[4+9*12]  = -D/s;

	    /* RpZZ */
	    MLR[5+0*12]  = -sp->k0*sp->k0*Deps;

	    /* ZpZZ */
	    MLR[11+6*12] = -sp->k0*sp->k0*Deps;

	    return;
	}

	/* No, we're inside the non-PML region. The following implements
	 * the interface equations (10)-(15) and (20)-(25) from the
	 * paper. (In wgms3d-1.0.0 and earlier versions, these were
	 * automatically generated using Maple, see the distribution
	 * archives of those versions.) */
	const double K = 1.0 + sp->c * rho;
	const complex<double> kkDeps = sp->k0 * sp->k0 * Deps;

	/* RpR */
	diffops_add_h1_entries(MLR+1, S*C*D);

	/* RpZ */
	diffops_add_h1_entries(MLR+2, S*S*D);

	/* RpRR */
	MLR[3+0*12] = -C*C*kkDeps; /* Coefficient for Rm */
	diffops_add_h1_entries(MLR+3, -S*D*(d*(4*C*C-1) + sp->c*C*C*C/K));
	diffops_add_h2_entries(MLR+3, C, S, -2.0*D*S*S*C);

	/* RpRZ */
	MLR[4+0*12] = -S*C*kkDeps; /* Coefficient for Rm */
	diffops_add_h1_entries(MLR+4, C*D*(d*(4*C*C-3) + sp->c*C*S*S/K));
	diffops_add_h2_entries(MLR+4, C, S, D*S*(2*C*C-1));

	/* RpZZ */
	MLR[5+0*12] = -S*S*kkDeps; /* Coefficient for Rm */
	diffops_add_h1_entries(MLR+5, S*D*(d*(4*C*C-1) - sp->c*C*S*S/K));
	diffops_add_h2_entries(MLR+5, C, S, 2.0*D*S*S*C);

	/* ZpR */
	diffops_add_h1_entries(MLR+7, -C*C*D);

	/* ZpZ */
	diffops_add_h1_entries(MLR+8, -S*C*D);

	/* ZpRR */
	MLR[9+6*12] = -C*C*kkDeps; /* Coefficient for Zm */
	diffops_add_h1_entries(MLR+9, C*D*(d*(4*C*C-3) + sp->c*C*C*C/K));
	diffops_add_h2_entries(MLR+9, C, S, 2.0*D*C*C*S);

	/* ZpRZ */
	MLR[10+6*12] = -S*C*kkDeps; /* Coefficient for Zm */
	diffops_add_h1_entries(MLR+10, S*D*(d*(4*C*C-1) + sp->c*C*C*C/K));
	diffops_add_h2_entries(MLR+10, C, S, -C*D*(2*C*C-1));

	/* ZpZZ */
	MLR[11+6*12] = -S*S*kkDeps; /* Coefficient for Zm */
	diffops_add_h1_entries(MLR+11, -C*D*(d*(4*C*C-3) - sp->c*C*S*S/K));
	diffops_add_h2_entries(MLR+11, C, S, -2.0*D*C*C*S);
    }

    void
    Diffops::do_matched_taylor_expansion (complex<double> *dstR,
					  complex<double> *dstZ,
					  int incd,
					  double rp,
					  double zp,
					  double dr,
					  double dz,
					  complex<double> epsp,
					  int &found_interfaces)
    {
	int N = 2*(NDO+1);
	std::list<nboundary> *interfaces
	    = mgp->find_intersections_with_line_segment(rp, zp, dr, dz);
	int M = 1;
	char trans = 'N';
	complex<double> alpha = 1.0, beta = 0.0;
	complex<double> *C = new complex<double>[N*N];
	complex<double> *D = new complex<double>[N*N];
	complex<double> *matrix = new complex<double>[N*N];

	std::memset(matrix, 0, sizeof(matrix[0])*N*N);
	int i;
	for(i = 0; i < N; i++) {
	    matrix[i+i*N] = 1.0;
	}

	if(debugwgms3d &&interfaces->size() > 0) {
	    std::cout << "We are at (" << rp << "," << zp << ")um, Delta=(" << dr << "," << dz << ")um" << std::endl;
	    std::cout << "number of interfaces = " << interfaces->size() << std::endl;
	}
    
	double lasta = 0.0;
	for(auto it = interfaces->begin(); it != interfaces->end();
	    lasta = it->a, epsp = it->epsr, ++it) {

	    if(debugwgms3d) {
		std::cout << " epsp is " << epsp << std::endl;
		std::cout << " Interface @a=" << it->a << " with c=" << it->c << "  " << it->epsl << " " << it->epsr << std::endl;
		std::cout << "   intersection point = " << rp+it->a*dr << " / " << zp+it->a*dz << std::endl;
	    }

	    if(debugmgp) {
		std::cout << " Interface at a=" << it->a << " with c=" << it->c << "  " << sqrt(it->epsl) << " " << sqrt(it->epsr)
			  << "   theta=" << it->theta << std::endl;
	    }

	    found_interfaces++;

	    if(it->a <= 1e-14 || it->a >= 1.0-1e-14) {
		std::cerr << "Grid problem: grid point right on dielectric interface at ("
			  << rp+it->a*dr << "," << zp+it->a*dz << ")" << std::endl;
		exit(1);
	    }

	    /* FIXME: sometimes the old lib2geom geometry handling reports
	     * multiple intersections even though there should be only
	     * one, I guess due to round-off errors. Try to detect these
	     * situations here. */
	    if(std::fabs(it->a - lasta) <= 1e-4 && it->epsl != epsp) {
		std::cout << "WARNING: Badly conditioned geometry." << std::endl;
		std::cout << "  Last intersection at (" << rp+it->a*dr << "," << zp+it->a*dz << ")." << std::endl;
		std::cout << "  We are at (x,y)=" << rp << "," << zp << "; (dr,dz)=" << dr << "," << dz << std::endl;
		std::cout << "  Complete intersection list:" << std::endl;
		for(auto it2 : *interfaces) {
		    it2.print();
		    /* I like C++11 */
		}
		continue;
	    }

	    if(it->epsl != epsp) {
		std::cout << "GEOMETRY ERROR: it->epsl=" << it->epsl << " doesn't match epsp=" << epsp << std::endl;
		std::cout << "  We are at (x,y)=" << rp << "," << zp << "; (dr,dz)=" << dr << "," << dz << std::endl;
		exit(1);
	    }

	    /* Set up matrix that expresses the unknown field and its
	     * derivatives (= vector f in the JLT paper) at point L (= in
	     * front of the interface) in terms of the vector f at the
	     * previous point. */
	    make_taylor_matrix(D, (it->a - lasta)*dr, (it->a - lasta)*dz);
	    /* C = D * matrix; */
	    GEMM(&trans, &trans, &N, &N, &N, &alpha, D, &N, matrix, &N, &beta, C, &N);

	    /* Set up matrix that expresses the vector f at point R (=
	     * behind the interface) in terms of the vector f at point L
	     * (= in front of the interface). In other words, apply the
	     * interface equations from Section III-B in the JLT paper. */
	    make_curv_interface_matrix(D, it->theta, it->c, it->epsl, it->epsr,
				       rp + it->a * dr, zp + it->a * dz);

	    /* matrix = D * C; */
	    GEMM(&trans, &trans, &N, &N, &N, &alpha, D, &N, C, &N, &beta, matrix, &N);
	}

	/* Now do the final 'Taylor step' to the desired end point, and
	 * store result for H_rho and H_z. */
	make_taylor_matrix(D, (1.0 - lasta)*dr, (1.0 - lasta)*dz);
	GEMM(&trans, &trans, &M, &N, &N, &alpha, D + 0,       &N, matrix, &N, &beta, dstR, &incd);
	GEMM(&trans, &trans, &M, &N, &N, &alpha, D + (NDO+1), &N, matrix, &N, &beta, dstZ, &incd);

	delete[] C;
	delete[] D;
	delete[] matrix;
	delete interfaces;
    }

    complex<double> *
    Diffops::get_standard_diffop (double n,
				  double e,
				  double s,
				  double w)
    {
	int k;

	std::memset(_stddiffop, 0, sizeof(_stddiffop));

	for(k = 0; k < 2; k++) {
#if NUM_GHOST_POINTS == 1
	    if(sp->use_five_point_standard) {
#include "standarddiffop_handmade.h"
	    } else {
#include "standarddiffop_2nd.h"
	    }
#endif

#if NUM_GHOST_POINTS == 2
	    exit(1);
#endif
	}

	return _stddiffop;
    }

/* get_diffops(): return previously calculated FD approximation to
 * differential operators at given point. Needed for
 * post-processing. */
    complex<double> *
    Diffops::get_diffops (const std::vector<double> &rs,
			  const std::vector<double> &zs,
			  int i,
			  int j)
    {
	complex<double> *d = diffops[(j << 16) + i];

	if(d == NULL) {
	    d = get_standard_diffop(zs[j+1]-zs[j], rs[i+1]-rs[i], zs[j]-zs[j-1], rs[i]-rs[i-1]);
	}

	return d;
    }


    Diffops::Diffops (MGP *waveguide_geometry,
		      std::shared_ptr<SimulationParameters> simulation_parameters)
    {
	TayA_M = 2*NDIRS;
	TayA_N = 2*NDO;
	Tay_nrhs = 2*(NDIRS+1);
	Tay_trans = 'N';
	Tay_lwork = -1;

	mgp = waveguide_geometry;
	sp = simulation_parameters;

	int info;
	complex<double> Tay_wwork0;

	/* Find out optimal work-array size for ZGELS */
	GELS(&Tay_trans, &TayA_M, &TayA_N, &Tay_nrhs, NULL, &TayA_M,
	     NULL, &TayA_M, &Tay_wwork0, &Tay_lwork, &info);
	if(info != 0) {
	    std::cerr << "ZGELS for LWORK=-1 failed with INFO = " << info << "." << std::endl;
	    exit(1);
	}

	Tay_lwork = (int)(Tay_wwork0.real());
	Tay_wwork = new complex<double>[Tay_lwork];
    }

    Diffops::~Diffops (void)
    {
	delete[] Tay_wwork;

	/* Free all stored diffops. (code quality TODO: don't store
	 * raw pointers; either smart pointers, or directly an
	 * array). */
	for_each(diffops.begin(), diffops.end(), delete_diffop);
    }

/* Calculate finite-difference approximations for differential
 * operators at point (rp,zp). The FD weights may be complex if we're
 * inside a PML region, since the interface conditions depend on the
 * complex stretching function s. However, even for points in the
 * non-PML regions, we return a complex<double> array in order to
 * minimize code duplication. */

/* Returns an array that must be freed with delete[] by the caller. */

    complex<double> *
    Diffops::calculate_diffop (double rp,
			       double zp,
			       complex<double> epsp,
			       const direction *dirs)
    {
	int k;
	int nrows, inc, info;
	complex<double> scale;
	complex<double> *M0;

	int found_interfaces = 0;

	/* Get Taylor expansions of fields to neighbouring
	 * mesh points: two field components, NDIRS
	 * neighbouring mesh points => 2*NDIRS expansions in
	 * terms of field and derivatives at P (=
	 * 2*(NDO+1)) */
    
	/* Tay is 2*NDIRS x 2*(NDO+1) */
    
	/* Tay will contain the Taylor expansions for the
	 * fields at the selected stencil points in terms of
	 * the field and its derivatives at the center point
	 * P. 
	 *
	 * The rows Tay(i,:) and Tay(NDIRS+i,:) (where 1 <= i <=
	 * NDIRS) contain the Taylor expansions for the i-th
	 * stencil point of H^r and H^z, respectively. In each
	 * row, the order of the coefficients is
	 * [ H^r, D_1 H^r, D_2 H^r, ... D_NDO H^r, 
	 *   H^z, D_1 H^z, D_2 H^z, ... D_NDO H^z ], so that
	 *
	 * H^r(i) = Tay(i,1)*H^r(P) + Tay(i,NDO+2)*H^z(P)
	 *          + \sum_{j=1}^{NDO} (Tay(i,j+1)*D_jH^r(i)
	 *                         + Tay(i,j+NDO+2)*D_jH^z(i))
	 * H^z(i) = Tay(i+NDO,1)*H^r(P) + Tay(i+NDO,NDO+2)*H^z(P)
	 *          + \sum_{j=1}^{NDO} (Tay(i+NDO,j+1)*D_jH^r(i)
	 *                         + Tay(i+NDO,j+NDO+2)*D_jH^z(i))
	 *
	 */

	complex<double> Tay[(2*NDIRS) * (2*(NDO+1))];
	std::memset(Tay, 0, sizeof(Tay));

	for(k = 0; k < NDIRS; k++) {
	    if(debugmgp) {
		std::cout << " direction #" << k << std::endl;
	    }
	    do_matched_taylor_expansion(Tay + k, Tay + k + NDIRS,
					(2*NDIRS), rp, zp, dirs[k].dr, dirs[k].dz, epsp,
					found_interfaces);
	}

	if(ENABLE_STANDARD_DIFFOPS && !found_interfaces) {
	    /* we really are in a homogeneous region. return NULL, let the
	     * main program use the explicit FD weights for that case */
	    return NULL;
	}
    
	/* Set up linear system of equations (system matrix is
	 * Tay*A, right-hand side is C) */
	complex<double> TayA[(2*NDIRS) * (2*NDO)];
	nrows = (2*NDIRS)*NDO;
	inc = 1;
	COPY(nrows, Tay + 1*(2*NDIRS), inc, TayA + 0*(2*NDIRS), inc);
	COPY(nrows, Tay + (NDO+2)*(2*NDIRS), inc, TayA + NDO*(2*NDIRS), inc);

	complex<double> C[(2*NDIRS) * ((2*NDIRS)+2)];
	std::memset(C, 0, sizeof(C));
	scale = -1.0;
	nrows = (2*NDIRS);
	inc = 1;
	COPY(nrows, Tay + 0*(2*NDIRS), inc, C + 0*(2*NDIRS), inc);
	SCAL(nrows, scale, C + 0*(2*NDIRS), inc);
	COPY(nrows, Tay + (NDO+1)*(2*NDIRS), inc, C + (NDIRS+1)*(2*NDIRS), inc);
	SCAL(nrows, scale, C + (NDIRS+1)*(2*NDIRS), inc);
	for(k = 0; k < NDIRS; k++) {
	    C[(k+1)*(2*NDIRS) + k] = 1.0;
	    C[(k+2+NDIRS)*(2*NDIRS) + k+NDIRS] = 1.0;
	}

	/* Compute least-squares solution */

	/* Solution expresses the (2*NDO) field derivatives
	 * in terms of fields at P and at those NDIRS
	 * neighbouring mesh points specified by the dirs[]
	 * array (this need not include all the stencil
	 * points!).

	 * TayA is 2*NDIRS x 2*NDO       = (M x N in xGELS)
	 * C    is 2*NDIRS x 2*(NDIRS+1) 

	 * (TayA^{-1}*C) is thus 2*NDO x 2*(NDIRS+1).

	 * The following calculates (TayA^{-1} * C) in a
	 * least-squares sense.

	 */
	GELS(&Tay_trans, &TayA_M, &TayA_N, &Tay_nrhs, TayA, &TayA_M, C, &TayA_M,
	     Tay_wwork, &Tay_lwork, &info);
	if(info != 0) {
	    std::cerr << "ZGELS failed with INFO = " << info << "." << std::endl;
	    exit(1);
	}

	/* Now sort into the M0 array.
       
	   while the upper result is 2*NDO x 2*(NDIRS+1),
	   (the column numbers refer to the directions in dirs[])
       
	   M0 is 2*NDO x 2*NSP
	   (the column numbers refer to the entire stencil)
       
	*/

	M0 = new complex<double>[(2*NDO) * (2*NSP)];
	std::memset(M0, 0, (2*NDO) * (2*NSP) * sizeof(*M0));

	nrows = 2*NDO;
	inc = 1;
	for(k = 0; k < NDIRS + 1; k++) {
	    /* (2*NDIRS) is the leading dimension of array C... */
	    COPY(nrows, C +          k *(2*NDIRS), inc,
		 M0 +      dirs_to_stencil_map[k] *(2*NDO), inc);
	    COPY(nrows, C + (NDIRS+1+k)*(2*NDIRS), inc,
		 M0 + (NSP+dirs_to_stencil_map[k])*(2*NDO), inc);
	}

	return M0;
    }

    void
    Diffops::store_diffops (complex<double> *M0,
			    int i,
			    int j)
    {
	/* We here merely store the pointer M0 (to an array of
	 * complex<double>). The arrays pointed to will be deleted when
	 * this Diffops instance is deleted. */
	diffops[(j << 16) + i] = M0;
    }

} // namespace wgms3d
