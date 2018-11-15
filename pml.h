
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

#ifndef _PML_H
#define _PML_H

#include <cassert>
#include <complex>

/* This class defines a PML.

   Corresponding to the real coordinate x, a complex path z(x) is
   defined according to

   z(x) = { x                                        for dir*(x-x0) < 0,
            x + j*strength*dir^m*(x-x0)^(m+1)/(m+1)  otherwise

   where x0 is the beginning of the PML, and 'dir' is its direction
   (+1 if it extends to positive x, and -1 if it extends towards
   negative x). The factor "dir^m" makes sure that the imaginary part
   of z(x) is always >= 0 (in order to damp (as opposed to amplify)
   outgoing waves).

   The stretching function is the derivative s = dz/dx:

   s(x) = { 1.0                                      for dir*(x-x0) < 0,
            1.0 + j*strength*dir^m*(x-x0)^m          otherwise

   We will also need the derivative of the stretching function:

   s'(x) = { 0.0                                     for dir*(x-x0) < 0,
             0.0 + j*strength*dir^m*m*(x-x0)^(m-1)   otherwise.

*/

class PML {
private:
    /* where the imaginary part of stretching function starts to
     * become non-zero: */
    double x0;
    /* whether PML extends towards positive (+1) or negative (-1)
     * direction. 0 means: PML disabled.*/
    int dir;

    double strength;
    int m;
    int pow_dir_m; /* = dir^m */

public:
    PML(void) {
	dir = 0;
	strength = 0.0;
    }

    void init (double x0,
	       int dir,
	       int m) {
	this->x0 = x0;
	this->dir = dir;
	this->m = m;
	if((m & 1) == 0) {
	    /* m even */
	    pow_dir_m = 1;
	} else {
	    /*m odd */
	    pow_dir_m = dir;
	}
    }

    /* Set optimal strength for PML. See Taflove/Hagness Eq. (7.61) */
    void set_optimal_strength (double h, /* grid-point distance */
			       double thickness, /* total thickness of PML */
			       double k0, /* free-space wave vector */
			       double multiplier) {
	assert(dir != 0); /* make sure init() has been called before. */
	this->strength = multiplier * (0.8 * (m+1)) / (k0 * h * pow(thickness,m));
    }

    /* check whether given point x is inside the PML. */
    bool is_inside (double x) const {
	if(dir == 0) {
	    /* no PML present. */
	    return false;
	} else {
	    return dir*(x-x0) > 0;
	}
    }

    /* return stretched coordinate z(x). */
    std::complex<double> get_stretched_x (double x) const {
	if(!is_inside(x)) {
	    return x;
	} else {
	    std::complex<double> z(x, strength * pow_dir_m * pow(x-x0, m+1) / (m+1));
	    return z;
	}
    }

    /* evaluate stretching function s(x). */
    std::complex<double> get_s (double x) const {
	if(!is_inside(x)) {
	    return 1.0;
	} else {
	    std::complex<double> s(1.0, strength * pow_dir_m * pow(x-x0, m));
	    return s;
	}
    }

    /* evaluate derivative of stretching function. we assume it is
     * differentiable at x=x0, i.e., m >= 2. */
    std::complex<double> get_sprime (double x) const {
	if(!is_inside(x)) {
	    return 0.0;
	} else {
	    std::complex<double> s(0.0, strength * pow_dir_m * m * pow(x-x0, m-1));
	    return s;
	}
    }
};

#endif
