
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

#ifndef _DIFFOPS_H
#define _DIFFOPS_H

#include <complex>
#include <unordered_map>

#include "wgms3d.h"
#include "stencil.h"
#include "mgp.h"

/* The diffop arrays (M0) contain the expressions for the derivatives
 * in terms of the surrounding mesh points.

 * M0(1,:) is expression for H^x_x
 * M0(2,:) is expression for H^x_y
 * M0(3,:) is expression for H^x_xx
 * M0(4,:) is expression for H^x_xy
 * M0(5,:) is expression for H^x_yy
 * etc.
 * M0(NDO+1,:) is expression for H^y_x
 * M0(NDO+2,:) is expression for H^y_y
 * etc.

 * The order of the coefficients in one row of M0 is:

 * M0(:,1:NSP) are for Hx values at stencil points (starting with P)
 * M0(:,NSP+1:(2*NSP)) are for Hy values at stencil points (starting with P)

 * M0 dimension is 2*NDO x 2*NSP

 */

namespace wgms3d {

    struct SimulationParameters;

    class Diffops
    {
    public:
	Diffops (MGP *waveguide_geometry,
		 std::shared_ptr<SimulationParameters> simulation_parameters);

	~Diffops (void);

	/* Calculate finite-difference approximations for differential
	 * operators at point (xp,yp). The FD weights may be complex if
	 * we're inside a PML region, since the interface conditions
	 * depend on the complex stretching function s. However, even for
	 * points in the non-PML regions, we return a complex<double>
	 * array in order to minimize code duplication. */
	/* Returns an array that must be freed with delete[] by the
	 * caller. */
	std::complex<double> * calculate_diffop (double xp,
						 double yp,
						 std::complex<double> epsp,
						 const direction *dirs);

	/* Returns an array that must not be freed -- it's only valid
	 * until the next call to this function: */
	std::complex<double> * get_standard_diffop (double n,
						    double e,
						    double s,
						    double w);

	void store_diffops (std::complex<double> *M0,
			    int i,
			    int j);

	/* get_diffops(): return previously calculated FD approximation to
	 * differential operators at given point. Needed for
	 * post-processing. */
	std::complex<double> * get_diffops (const std::vector<double> &rs,
					    const std::vector<double> &zs,
					    int i,
					    int j);

	int get_num_stored_diffops (void) {
	    return diffops.size();
	}

    private:
	int TayA_M, TayA_N, Tay_nrhs;
	char Tay_trans;
	int Tay_lwork;
	std::complex<double> *Tay_wwork;

	std::unordered_map<int, std::complex<double>*> diffops;

	std::complex<double> _stddiffop[(2*NDO) * (2*NSP)];

	MGP *mgp;
	std::shared_ptr<SimulationParameters> sp;

	void make_curv_interface_matrix (
	    std::complex<double> *MLR,
	    double theta, /* angle of interface at intersection point */
	    double d, /* curvature of interface */
	    std::complex<double> m, /* n^2 before interface */
	    std::complex<double> p, /* n^2 behind interface */
	    double rho, double z); /* non-stretched (real) coordinates of intersection point*/

	/* This function returns the field values H^{rho|z}(rp+dr,zp+dz)
	   in terms of the values of the fields H^{rho|z} and their
	   derivatives at (rp,zp). */
	void do_matched_taylor_expansion (
	    std::complex<double> *dstR,
	    std::complex<double> *dstZ,
	    int incd,
	    double rp,
	    double zp,
	    double dr,
	    double dz,
	    std::complex<double> epsp,
	    int &found_interfaces);

    };

} // namespace wgms3d

#endif // _DIFFOPS_H
