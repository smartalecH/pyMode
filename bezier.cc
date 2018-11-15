
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>

 
    This file contains some class definitions and code from lib2geom
    snapshot 20101112, written and Copyright (C) by Marco Cecchetti
    <mrcekets at gmail.com>, Michael G. Sloan <mgsloan@gmail.com>, and
    others, and originally released under LGPL v2.1.

    Polynomial Root-Finder Copyright (c) 2003, by Per Vognsen. Adapted
    from http://www.flipcode.com/archives/Polynomial_Root-Finder.shtml
    on 20110204. Original license: "This software is free for any
    use."


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
#include <iterator>
#include <vector>
#include <algorithm>
#include <cassert>

#include "bezier.h"

static Point
lerp (double const t,
      Point const &a,
      Point const &b)
{
    return (a * (1 - t) + b * t);
}

/*
 * Compute the hodograph of the bezier curve B and return it in D
 *
 * See
 * http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html:
 * "the derivative of C(u) is a Bézier curve of degree n - 1 defined
 * by n control points n(P1 - P0), n(P2 - P1), n(P3 - P2), ..., n(Pn -
 * Pn-1)"
 */
void
bezier_Derivative (std::vector<Point> & D,
		   std::vector<Point> const& B)
{
    D.clear();
    size_t sz = B.size();
    if (sz == 0) return;
    if (sz == 1)
    {
        D.resize(1, Point(0,0));
        return;
    }
    size_t n = sz-1;
    D.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        D.push_back((B[i+1] - B[i])*n);
    }
}

static Coord
dot (Point const &a,
     Point const &b)
{
    return a[0] * b[0] + a[1] * b[1];
}

/*
 *  Compute the portion of the Bezier curve "B" wrt the interval [t,1]
 *  (de Casteljau's algorithm)
 */
static void
right_portion (Coord t,
	       std::vector<Point> & B)
{
    size_t n = B.size();
    for (size_t i = 1; i < n; ++i)
    {
        for (size_t j = 0; j < n-i; ++j)
        {
            B[j] = lerp(t, B[j], B[j+1]);
        }
    }
}

Point
bezier_valueAt (std::vector<Point> const &B,
		double t)
{
    std::vector<Point> C = B;
    right_portion(t, C);
    return C[0];
}

template<typename T>
int sign (T x)
{
    if (x > 0)
	return 1;
    else if (x < 0)
	return -1;
    else
	return 0;
}

static double
evaluate (const std::vector<double>& coeff,
	  double x)
{
    assert(!coeff.empty());
	
    double value = 0, y = 1;
    for (unsigned int i = 0; i < coeff.size(); i++) {
	value += coeff[i] * y;
	y *= x;
    }
    return value;
}

static double
bisect (const std::vector<double>& coeff,
	double low,
	double high)
{
    /* 100 iterations means that our initial search-interval size of
     * about unity will be reduced to 1/2^100 here. That accuracy is
     * surely enough for our purposes. */
    const unsigned int max_iterations = 100;
    double mid;

    for(unsigned int i = 0; i < max_iterations; i++) {
	mid = 0.5 * (low + high);
	if(sign(evaluate(coeff, low)) != sign(evaluate(coeff, mid)))
	    high = mid;
	else
	    low = mid;
    }

    return mid;
}
	
static std::vector<double>
differentiate (const std::vector<double>& coeff)
{
    if (coeff.size() < 2)
	return std::vector<double>();

    std::vector<double> deriv;
    deriv.resize(coeff.size()-1);
    for (unsigned int i = 0; i < deriv.size(); i++)
	deriv[i] = (i+1) * coeff[i+1];
    return deriv;
}

/*
  Return all real roots of polynomial in interval [0,1] (plus possibly
  some more roots outside this interval).
 */

static std::vector<double>
find_roots (const std::vector<double>& coeff)
{
    assert(coeff.size() >= 2);

#if 0
    std::cout << "Polynomial: ";
    std::copy(coeff.begin(), coeff.end(),
	      std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
#endif

    std::vector<double> roots;
    if(coeff.size() == 2) {
	if(coeff[1] == 0.0) {
	    if(coeff[0] == 0.0) {
		std::cerr << "ERROR: find_roots(): infinite number of zeros." << std::endl;
		exit(1);
	    } else {
		/* no root */
		return roots;
	    }
	} else {
	    roots.push_back(-coeff[0]/coeff[1]);
	}
    } else {
	std::vector<double> deriv_roots(find_roots(differentiate(coeff)));

	if (deriv_roots.empty())
	    deriv_roots.push_back(0.0);

	/* We're ultimately interested in roots inside [0,1] only, so
	 * limit the search intervals right here. */
	std::vector<double>::iterator it = deriv_roots.begin();
	while(it < deriv_roots.end()) {
	    if(*it < -0.5 || *it > +1.5) {
		deriv_roots.erase(it);
		it = deriv_roots.begin();
	    } else {
		it++;
	    }
	}
	deriv_roots.push_back(-0.6);
	deriv_roots.push_back(+1.6);

	std::sort(deriv_roots.begin(), deriv_roots.end());
	deriv_roots.erase(std::unique(deriv_roots.begin(), deriv_roots.end()),
			  deriv_roots.end());

	for (unsigned int i = 0; i < deriv_roots.size()-1; i++) {
	    if(   sign(evaluate(coeff, deriv_roots[i])) 
	       != sign(evaluate(coeff, deriv_roots[i+1]))) {
		roots.push_back(bisect(coeff, deriv_roots[i], deriv_roots[i+1]));
	    }
	}

	std::sort(roots.begin(), roots.end());
	roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    }

#if 0
    std::cout << "Roots found: ";
    std::copy(roots.begin(), roots.end(),
	      std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
#endif

    return roots;
}

void
bezier_find_line_intersections (std::vector< std::pair<double, double> >& xs,
				std::vector<Point> const& bezier,
				double px,
				double py,
				double dx,
				double dy)
{
    /* Get normal form of line: x * n = d */
    Point n(dy,-dx);
    double d = px*n[0] + py*n[1];

    std::vector<double> coeff;

    if(bezier.size() == 5) {
	/* Construct polynomial that describes intersection of Bézier
	 * curve and line. This is for the special case that A, B, C are
	 * the three control points specified by the user and that two
	 * more control points are interpolated between these three to
	 * make the local curvature at the ends of the curve zero (see
	 * mgp.cc). */
	Point A = bezier[0];
	Point B = bezier[2];
	Point C = bezier[4];

	Point G = (B-A)*2;
	Point F = (A+C-B*2)*2;
	Point E = B*2-A-C;

	coeff.push_back(dot(A,n)-d);
	coeff.push_back(dot(G,n));
	coeff.push_back(0.0);
	coeff.push_back(dot(F,n));
	coeff.push_back(dot(E,n));
    } else if(bezier.size() == 3) {
	/* Construct polynomial that describes intersection of Bézier
	 * curve and line. This is for the special case that A, B, C
	 * are the three control points specified by the user and no
	 * other points are interpolated. */
	Point A = bezier[0];
	Point B = bezier[1];
	Point C = bezier[2];

	coeff.push_back(dot(A,n)-d);
	coeff.push_back(dot(B*2 - A*2,n));
	coeff.push_back(dot(A - B*2 + C,n));
    } else {
	std::cerr << "PROGRAMMING ERROR in Bézier handling." << std::endl;
	exit(1);
    }

    /* Remove leading-zero coefficients */
    while(coeff.size() > 0 && coeff[coeff.size()-1] == 0.0) {
	coeff.resize(coeff.size()-1);
	//std::cout << "NOTE: reducing degree of Bézier-intersection polynomial." << std::endl;
    }
    if(coeff.size() < 2) {
	std::cerr << "GEOMETRY ERROR: degenerate Bézier-intersection polynomial." << std::endl;
	exit(1);
    }

    /* Find intersection points */
    std::vector<double> roots = find_roots(coeff);
    for(unsigned int i = 0; i < roots.size(); i++) {
	if(roots[i] >= 0.0 && roots[i] <= 1.0) {
	    Point Bu = bezier_valueAt(bezier, roots[i]);

	    /* Calculate value of t (fraction of line segment where
	     * intersection happens) */
	    double t = ((Bu[0] - px) * dx + (Bu[1] - py) * dy) / (dx*dx + dy*dy);

	    /* If t is in [0,1], we have an intersection; add to
	     * output list */
	    if(t >= 0.0 && t <= 1.0) {
		std::pair<double, double> ci;		
		ci.first = roots[i];
		ci.second = t;
		xs.push_back(ci);
	    }
	}
    }
}
