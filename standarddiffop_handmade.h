
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

	/* \partial rho */
	double tmp = w*e * (w+e);
	_stddiffop[0+NDO*k + (3+NSP*k)*(2*NDO)] = w*w / tmp;
	_stddiffop[0+NDO*k + (7+NSP*k)*(2*NDO)] = -e*e / tmp;
	_stddiffop[0+NDO*k + (0+NSP*k)*(2*NDO)] = (e*e - w*w) / tmp;

	/* \partial z */
	tmp = s*n * (s+n);
	_stddiffop[1+NDO*k + (1+NSP*k)*(2*NDO)] = s*s / tmp;
	_stddiffop[1+NDO*k + (5+NSP*k)*(2*NDO)] = -n*n / tmp;
	_stddiffop[1+NDO*k + (0+NSP*k)*(2*NDO)] = (n*n - s*s) / tmp;

	/* \partial rho \partial rho */
	tmp = 2.0 / (w*e * (w+e));
	_stddiffop[2+NDO*k + (3+NSP*k)*(2*NDO)] = w * tmp;
	_stddiffop[2+NDO*k + (7+NSP*k)*(2*NDO)] = e * tmp;
	_stddiffop[2+NDO*k + (0+NSP*k)*(2*NDO)] = -(w+e) * tmp;

	/* \partial z \partial z */
	tmp = 2.0 / (s*n * (s+n));
	_stddiffop[4+NDO*k + (1+NSP*k)*(2*NDO)] = s * tmp;
	_stddiffop[4+NDO*k + (5+NSP*k)*(2*NDO)] = n * tmp;
	_stddiffop[4+NDO*k + (0+NSP*k)*(2*NDO)] = -(n+s) * tmp;

	/* \partial rho \partial z */
	/* Note that the mixed partial derivative is only used for
	 * calculating the E field from the H field, it's not used in
	 * the actual mode solving process. */
	tmp = 1.0 / (n*e + w*s + s*e + w*n);
	_stddiffop[3+NDO*k + (2+NSP*k)*(2*NDO)] = tmp;
	_stddiffop[3+NDO*k + (4+NSP*k)*(2*NDO)] = -tmp;
	_stddiffop[3+NDO*k + (6+NSP*k)*(2*NDO)] = tmp;
	_stddiffop[3+NDO*k + (8+NSP*k)*(2*NDO)] = -tmp;
