
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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cassert>
using std::ifstream;

#include "mgp.h"
#include "bezier.h"

namespace {

    double norm_angle (double a)
    {
	while(a > M_PI) {
	    a -= 2*M_PI;
	}
	while(a < -M_PI) {
	    a += 2*M_PI;
	}
	return a;
    }

} // anonymous namespace

namespace wgms3d {

    std::list<nboundary> *
    MGP::find_intersections_with_line_segment (double px,
					       double py,
					       double dx,
					       double dy) const
    {
	double theta;
	std::complex<double> nm, np;
	std::list<nboundary> *nbs = new std::list<nboundary>;

	std::vector<Point> ptsLine;
	ptsLine.push_back(Point(px, py));
	ptsLine.push_back(Point((px+dx), (py+dy)));
	Rect rayboundingbox = Rect::from_range(ptsLine.begin(), ptsLine.end());

	for(auto lit = lin_interfaces.begin(); lit != lin_interfaces.end(); ++lit) {
	
	    double D0 = dx*(lit->y1-lit->y2) - dy*(lit->x1-lit->x2);
	    if(D0 != 0.0) {
		double Da = (lit->x1-px)*(lit->y1-lit->y2)
		    - (lit->y1-py)*(lit->x1-lit->x2);
		double aa = Da / D0;
		if(aa >= 0.0 && aa <= 1.0) {
		    double Db = dx*(lit->y1-py) - dy*(lit->x1-px);
		    double bb = Db / D0;

		    if(bb >= 0.0 && bb <= 1.0) {
			theta = atan2(lit->y2-lit->y1, lit->x2-lit->x1) - M_PI/2.0;
			theta = norm_angle(theta);
		  
			/* Which way is the interface crossed? */
			if(D0 < 0) {
			    nm = lit->nl;
			    np = lit->nr;
			} else {
			    nm = lit->nr;
			    np = lit->nl;
			}
  
			nboundary ni = { aa, theta, 0.0, nm*nm, np*np };
		    
			nbs->push_back(ni);
		    }
		}
	    }
	}

	for(auto bit = bez_interfaces.begin(); bit != bez_interfaces.end(); ++bit) {

	    /* quick & rough bounding-box check */
	    if(!rayboundingbox.intersects(bit->boundingbox)) {
		continue;
	    }

	    std::vector< std::pair<double, double> > xs;
	    bezier_find_line_intersections(xs, bit->points, px, py, dx, dy);

	    for(size_t i = 0; i < xs.size(); ++i) {
		Point tangente = bezier_valueAt(bit->deriv1, xs[i].first);
		Point deriv2 = bezier_valueAt(bit->deriv2, xs[i].first);
		double c = (tangente[0]*deriv2[1] - deriv2[0]*tangente[1]) /
		    pow(pow(tangente[0],2.0) + pow(tangente[1],2.0), 1.5);
		// curvature left is positive
		// curvature right is negative

		theta = atan2(tangente[1],tangente[0]) + M_PI/2.0;
		if(c > 0) {
		    theta = theta + M_PI;
		} else {
		    c = -c;
		}
		theta = norm_angle(theta);

#if 0
		std::cout << "Intersection #" << i << " at " << xs[i].first << " ; " << xs[i].second << std::endl;
		std::cout << "  Curvature is " << c << std::endl;
		std::cout << "  Angle is " << theta << std::endl;
#endif

		/* find out epsL and epsR */
		double nx = dy;
		double ny = -dx;
		if((tangente[0] * nx + tangente[1] * ny) < 0) {
		    nm = bit->nl;
		    np = bit->nr;
		} else {
		    nm = bit->nr;
		    np = bit->nl;
		}
  
		nboundary ni = { xs[i].second, theta, c, nm*nm, np*np };
		nbs->push_back(ni);
	    }
	}

	for(auto eit = ell_interfaces.begin(); eit != ell_interfaces.end(); ++eit) {
	    double aa = pow(dx*eit->ry,2.0) + pow(dy*eit->rx,2.0);
	    double xp = px - eit->x;
	    double yp = py - eit->y;
	    double p = (2.0*xp*dx*pow(eit->ry,2.0) + 2.0*yp*dy*pow(eit->rx,2.0)) / aa;
	    double q = (pow(xp*eit->ry,2.0) + pow(yp*eit->rx,2.0) - pow(eit->rx*eit->ry,2.0)) / aa;
	    double D = p*p/4.0 - q;

	    if(D < 0) {
		continue;
	    }
	    double sqrtD = sqrt(D);
	    double a = -p/2.0 + sqrtD;
	    if(a < 0.0 || a > 1.0) {
		a = -p/2.0 - sqrtD;
	    }
	    if(a < 0.0 || a > 1.0) {
		continue;
	    }

	    double x0_minus_c_over_r2 = (xp + a*dx - eit->x) / pow(eit->rx, 2.0);
	    double y0_minus_d_over_p2 = (yp + a*dy - eit->y) / pow(eit->ry, 2.0);
	    theta = atan2(y0_minus_d_over_p2, x0_minus_c_over_r2);
	    double Rloc = pow(eit->rx * eit->ry, 2.0)
		* pow(pow(x0_minus_c_over_r2, 2.0)
		      + pow(y0_minus_d_over_p2, 2.0), 1.5);

	    if(q > 0) {
		nm = eit->nout;
		np = eit->nin;
	    } else if(q < 0) {
		nm = eit->nin;
		np = eit->nout;
	    } else {
		std::cerr << std::endl;
		std::cerr << "GEOMETRY ERROR: Elliptical interface hit." << std::endl;
		exit(1);
	    }

	    nboundary ni = { a, theta, 1.0/Rloc, nm*nm, np*np };
	    nbs->push_back(ni);
	}

	nbs->sort();

	return nbs;
    }

    std::unique_ptr< std::complex<double>[] >
    MGP::get_epsis_on_grid (const std::vector<double> &x,
			    const std::vector<double> &y) const
    {
	int ldepsis = x.size();
	std::unique_ptr< std::complex<double>[] > epsis( new std::complex<double>[ldepsis * y.size()] );

	/* Vertical scanning from top using default index n_def */
	int i, j;
	std::complex<double> epsi;
	double y0;
	std::list<nboundary> *ifs;

	for(i = 0; i < int(x.size()); i++) {
	    epsi = n_def * n_def;
	    /* starting y0 for scanning: */
	    y0 = y[y.size()-1] + (y[y.size()-1] - y[y.size()-2]);
	    if(scany0_set) {
		y0 = scany0;
	    }
	    for(j = y.size()-1; j >= 0; j--) {
		ifs = find_intersections_with_line_segment(x[i], y0, 0, y[j]-y0);
		if(ifs->size() > 0) {
		    const nboundary &iff = ifs->back();
		    epsi = iff.epsr;
		    delete ifs;
		}
		epsis[ldepsis*j + i] = epsi;
		y0 = y[j];
	    }
	}

	return epsis;
    }

    MGP::MGP (const char *fn)
    {
	n_def = 0.0;
	scany0_set = false;
	scany0 = 0.0;
        has_complex_indices = false;

	core_x1 = 0.0; core_y1 = 0.0; core_x2 = 0.0; core_y2 = 0.0;

	ifstream ifs (fn, ifstream::in);

	if(ifs.fail()) {
	    std::cerr << std::endl;
	    std::cerr << "Couldn't open geometry file '" << fn << "'" << std::endl;
	    exit(1);
	}

	while(!ifs.eof()) {
	    char a;

	    ifs >> a;	
	    if(!ifs.good()) {
		std::cerr << std::endl;
		std::cerr << "Error reading " << fn << std::endl;
		exit(1);
	    }

	    switch(a) {
	    case '#':
		/* Comment line */
		ifs.ignore(31337, '\n');
		break;

	    case 'x':
		/* Exit */
		goto ende;
		break;

	    case 'n':
		/* Background refractive index */
		ifs >> n_def;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (n)." << std::endl;
		    exit(1);
		}
		updateHasComplexIndices(n_def);
		break;

	    case 'y':
		/* Set scany0 */
		ifs >> scany0;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (y)." << std::endl;
		    exit(1);
		}
		scany0_set = true;
		break;

	    case 'C':
		/* Define core-region rectangle */
		ifs >> core_x1 >> core_y1 >> core_x2 >> core_y2;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (C)." << std::endl;
		    exit(1);
		}
		break;

	    case 'l': {
		/* New interface, linear */
		lin_interface nil;
		ifs >> nil.nl >> nil.nr >> nil.x1 >> nil.y1 >> nil.x2 >> nil.y2;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (l)." << std::endl;
		    exit(1);
		}
		updateHasComplexIndices(nil.nl);
		updateHasComplexIndices(nil.nr);
		lin_interfaces.push_back(nil);
		break;
	    }

	    case 'b': {
		/* New interface, three-point Bézier curve */
		bez_interface nib;
		double x1, y1, x2, y2, x3, y3;
		ifs >> nib.nl >> nib.nr >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (b)." << std::endl;
		    exit(1);
		}
		updateHasComplexIndices(nib.nl);
		updateHasComplexIndices(nib.nr);

		/* Even though the user specifies three Bézier control
		   points only, we here interpolate two more nodes between
		   points 1+2 and 2+3 to make sure the local curvature is
		   zero at the ends of the curve. */
		nib.points.push_back(Point(x1, y1));
		nib.points.push_back(Point((x1+x2)/2.0, (y1+y2)/2.0));
		nib.points.push_back(Point(x2, y2));
		nib.points.push_back(Point((x2+x3)/2.0, (y2+y3)/2.0));
		nib.points.push_back(Point(x3, y3));
		bezier_Derivative(nib.deriv1, nib.points);
		bezier_Derivative(nib.deriv2, nib.deriv1);
		nib.boundingbox = Rect::from_range(nib.points.begin(), nib.points.end());
		bez_interfaces.push_back(nib);
		break;
	    }

	    case 'c': {
		/* New interface, circular */
		ell_interface eic;
		ifs >> eic.nout >> eic.nin >> eic.x >> eic.y >> eic.rx;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (c)." << std::endl;
		    exit(1);
		}
		updateHasComplexIndices(eic.nout);
		updateHasComplexIndices(eic.nin);
		eic.ry = eic.rx;
		ell_interfaces.push_back(eic);
		break;
	    }
	    
	    case 'e': {
		/* New interface, elliptic */
		ell_interface eic;
		ifs >> eic.nout >> eic.nin >> eic.x >> eic.y >> eic.rx >> eic.ry;
		if(!ifs.good()) {
		    std::cerr << std::endl;
		    std::cerr << "Error parsing MGP file (e)." << std::endl;
		    exit(1);
		}
		updateHasComplexIndices(eic.nout);
		updateHasComplexIndices(eic.nin);
		ell_interfaces.push_back(eic);
		std::cerr << std::endl;
		std::cerr << "WARNING: Implementation of elliptic interfaces not yet verified." << std::endl;
		std::cerr << "         Use at your own risk." << std::endl;
		std::cerr << std::endl;
		break;
	    }

	    default:
		std::cerr << std::endl;
		std::cerr << "Error parsing MGP file, unknown code '" << a << "'" << std::endl;
		exit(1);
		break;
	    }
	}

      ende:
	ifs.close();
    }

    void MGP::updateHasComplexIndices (std::complex<double> new_index)
    {
	if(new_index.imag() != 0.0) {
	    has_complex_indices = true;
	}
    }

    bool MGP::hasComplexIndices () const
    {
	return has_complex_indices;
    }

} // namespace wgms3d
