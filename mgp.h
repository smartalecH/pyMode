
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

#ifndef _MGP_H
#define _MGP_H

#include <list>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>

#include "geometry.h"

namespace wgms3d {

    struct nboundary {
	double a;
	double theta;
	double c; /* curvature */
	std::complex<double> epsl, epsr;

	bool operator< (nboundary &b) {
	    return a < b.a;
	}

	void print (void) {
	    std::cout << "  Interface @a=" << a << " with c=" << c
		      << ", " << std::sqrt(epsl) << "->" << std::sqrt(epsr)
		      << ", theta=" << theta << std::endl;
	}
    };

    struct lin_interface {
	std::complex<double> nl, nr;
	double x1, y1, x2, y2;
    };

    struct bez_interface {
	std::complex<double> nl, nr;
	std::vector<Point> points;
	std::vector<Point> deriv1;
	std::vector<Point> deriv2;
	Rect boundingbox;
    };

    struct ell_interface {
	std::complex<double> nout, nin;
	double x, y, rx, ry;
    };

    class MGP
    {
    public:

	MGP (const char *mgp_filename);

	bool is_core_layer_defined (void) const
	{
	    return core_x1 != core_x2;
	}
    
	bool is_in_core (double x,
			 double y) const
	{
	    if(core_x1 == core_x2) {
		return 0;
	    }
	
	    return x >= core_x1 && x <= core_x2 && y >= core_y1 && y <= core_y2;
	}

	std::list<nboundary> *
	find_intersections_with_line_segment (double px,
					      double py,
					      double dx,
					      double dy) const;

	/// Return a (Fortran-style) array representing the relative
	/// permittivity (epsilon) on the specified grid. Leading
	/// dimension will be rho.size().
	std::unique_ptr< std::complex<double>[] >
	get_epsis_on_grid (
	    const std::vector<double> &rho,
	    const std::vector<double> &z
	    ) const;

	/// Returns true if at least one of the materials used in the
	/// structure has a complex (not purely real) effective index.
	bool hasComplexIndices () const;

    private:

	std::complex<double> n_def;
	bool scany0_set;
	double scany0;

	double core_x1, core_y1, core_x2, core_y2;

	std::list<lin_interface> lin_interfaces;
	std::list<bez_interface> bez_interfaces;
	std::list<ell_interface> ell_interfaces;

	bool has_complex_indices;

	void updateHasComplexIndices (std::complex<double> new_index);

    };

} // namespace wgms3d

#endif /* _MGP_H */
