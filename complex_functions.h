
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

#ifndef WGMS3D_COMPLEX_FUNCTIONS_H
#define WGMS3D_COMPLEX_FUNCTIONS_H

/** \file
 *
 * Some helper functions for dealing with complex numbers.
 */

#include <complex>

namespace wgms3d {

    /// Convert complex number to complex number (= do nothing),
    /// required for complex calculation mode.
    inline const std::complex<double> & maybeConvertComplexToReal (
	const std::complex<double> &x,
	std::complex<double> /* tag */)
    {
	return x;
    }

    /// Convert complex number to real number (= discarding the
    /// imaginary part), required for real calculation mode.
    inline double maybeConvertComplexToReal (
	const std::complex<double> &x,
	double /* tag */)
    {
	return x.real();
    }

} // namespace wgms3d
    
#endif // WGMS3D_COMPLEX_FUNCTIONS_H
