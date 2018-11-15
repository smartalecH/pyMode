
/*
    wgms3d - a full-vectorial finite-difference mode solver.

    Copyright (C) 2005-2012  Michael Krause <m.krause@tu-harburg.de>

 
    This file contains some class definitions and code from lib2geom
    snapshot 20101112, written and Copyright (C) by Marco Cecchetti
    <mrcekets at gmail.com>, Michael G. Sloan <mgsloan@gmail.com>, and
    others, and originally released under LGPL v2.1.


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

#ifndef __GEOMETRY_H
#define __GEOMETRY_H

#include <vector>

typedef double Coord;

class Point {
  private:
    Coord _pt[2];

  public:
    Point (Coord x,
	   Coord y) {
        _pt[0] = x; _pt[1] = y;
    }

    Coord operator[] (unsigned i) const {
	return _pt[i];
    }

    Point &operator+= (Point const &o) {
        for ( unsigned i = 0 ; i < 2 ; ++i ) {
            _pt[i] += o._pt[i];
        }
        return *this;
    }

    const Point operator+ (const Point &other) const {
	Point result = *this;
	result += other;
	return result;
    }

    Point &operator-= (Point const &o) {
        for ( unsigned i = 0 ; i < 2 ; ++i ) {
            _pt[i] -= o._pt[i];
        }
        return *this;
    }

    const Point operator- (const Point &other) const {
	Point result = *this;
	result -= other;     
	return result;       
    }

    Point operator-() const {
        return Point(-_pt[0], -_pt[1]);
    }

    Point &operator*= (Coord s) {
        for ( unsigned i = 0 ; i < 2 ; ++i ) _pt[i] *= s;
        return *this;
    }

    const Point operator* (const Coord &s) const {
	Point result = *this;
	result *= s;
	return result;       
    }
};

inline std::ostream &operator<< (std::ostream &out_file, const Point &in_pnt) {
    out_file << "X: " << in_pnt[0] << "  Y: " << in_pnt[1];
    return out_file;
}

class Interval {
  private:
    Coord _b[2];

  public:
    explicit Interval() { _b[0] = 0;  _b[1] = 0; }

    Interval(Coord u, Coord v) {
        if (u <= v) {
            _b[0] = u; _b[1] = v;
        } else {
            _b[0] = v; _b[1] = u;
        }
    }

    Coord operator[] (unsigned i) const { return _b[i]; }

    void expandTo(Coord val) {
       if(val < _b[0]) _b[0] = val;
       if(val > _b[1]) _b[1] = val;  //no else, as we want to handle NaN
    }

    /** @brief Check whether the interval includes this number. */
    bool contains(Coord val) const { return _b[0] <= val && val <= _b[1]; }

    /** @brief Check whether the interval includes the given interval. */
    bool contains(Interval const &val) const { return _b[0] <= val._b[0] && val._b[1] <= _b[1]; }

    bool intersects(Interval const &val) const {
        return contains(val._b[0]) || contains(val._b[1]) || val.contains(*this);
    }
};

class Rect {
private:
    Interval f[2];

public:
    Rect() {
	f[0] = f[1] = Interval();
    }

    Rect (Point const & a,
	  Point const & b) {
        f[0] = Interval(a[0], b[0]);
        f[1] = Interval(a[1], b[1]);
    }

    template <typename InputIterator>
    static Rect from_range(InputIterator start, InputIterator end) {
        assert(start != end);
        Point p1 = *start++;
        Rect result(p1, p1);
        for (; start != end; ++start) {
            result.expandTo(*start);
        }
        return result;
    }

    void expandTo(Point p) { 
        f[0].expandTo(p[0]);
	f[1].expandTo(p[1]);
    }

    bool intersects(Rect const &r) const { 
        return f[0].intersects(r[0]) && f[1].intersects(r[1]);
    }

    Interval const & operator[] (unsigned i) const {
	return f[i];
    }
};

#endif // __GEOMETRY_H
