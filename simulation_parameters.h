
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

/** \file
 *
 * Internal data structure for storing the simulation parameters.
 *
 */

#ifndef _WGMS3D_SIMULATION_PARAMETERS_H
#define _WGMS3D_SIMULATION_PARAMETERS_H

#include <vector>
#include <memory>

#include "wgms3d.h"
#include "mgp.h"
#include "pml.h"

namespace wgms3d
{

    /// Configuration of the mode solver. Will also be passed to
    /// instances of class Mode.
    struct SimulationParameters
    {
    public:
	FDMode fd_mode;
	bool use_five_point_standard;
	double k0; /* free-space wave number */
	double c; /* Waveguide curvature c = 1/R */
	PML pml[4]; /* pml_north, pml_east, pml_south, pml_west */

	/* Order: left, right, top, bottom. 0 = Electric wall, 1 =
	   Magnetic wall. */
	int bconds[4];

	/* 0 = grid point lies right on wall (default), 1 = two grid
	   points lie symmetrically around wall (for pseudo-2D mode,
	   currently non-functional) */
	int bcondsym[4];

	/* nir = 'n'umber of 'i'nner 'r'ho-grid points = all points that are
	 * not ghost points. */
	/// TODO get rid of this
	int nir, niz;

	/* number of values stored for a single field component: (=nir*niz) */
	/// TODO get rid of this
	int fcsize;

	/// Grid including ghost points
	/// TODO rename
	std::vector< double > _rhos, _zs;

	/// Grid (_rhos and _zs) with complex stretching
	std::vector< std::complex<double> > stretched_rhos, stretched_zs;

	/// Waveguide geometry
	std::unique_ptr< MGP > mgp;

	SimulationParameters (void) {
	    fd_mode = FDMode::FullVectorial;
	    use_five_point_standard = false;
	    k0 = 2 * M_PI / 1.55e-6;
	    c = 0;

	    nir = -1; niz = -1; fcsize = -1;
	    std::memset(bconds, 0, sizeof(bconds));
	    std::memset(bcondsym, 0, sizeof(bcondsym));
	}

	/// Returns an array 'retain_list' of size 2*nir*niz, one
	/// entry for each discretization point. If retain_list[i] <
	/// 0, then point i is known to be zero in advance (Dirichlet
	/// BC, or SV/scalar approximation).  Otherwise,
	/// retain_list[i] contains a non-negative integer. The
	/// non-negative entries of retain_list[] are consecutive
	/// integers starting at zero. TODO: required for non-root
	/// process?
	const std::vector<int> & getRetainList() const;

	/// 'number_of_unknowns' is the number of non-negative entries
	/// in retain_list.
	unsigned getNumUnknowns() const;

    private:
	mutable std::vector<int> retain_list;
	mutable unsigned number_of_unknowns;
	mutable std::mutex mutables;

	void ensureRetainList () const;

    };


} // namespace wgms3d

#endif // _WGMS3D_SIMULATION_PARAMETERS_H
